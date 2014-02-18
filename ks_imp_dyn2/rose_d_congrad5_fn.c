/******* d_congrad5_fn_tmp.c - conjugate gradient for SU3/fermions ****/
/* TEST VERSION  4/18/03, TIMING FOR GLOBAL REDUCTIONS */
/* REDUCE NUMBER OF GLOBAL SUMS */
/* MIMD version 6 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  
   This code combines d_congrad5_fn.c and d_congrad5_fn_tmp.c
   With CONGRAD_TMP_VECTORS defined, allocates temporary CG vectors in
   field-major order and uses them instead of the site-major
   temporaries.  They may be eliminated from the site structure in the
   future.
   Calls dslash_fn or dslash_fn_on_temp depending accordingly. */
/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.
	This is different than our old definition of the stopping
	criterion.  To convert an old stopping residual to the new
	one, multiply the old one by sqrt( (2/3)/(8+2*m) )
        This is because the source is obtained from
        a random vector with average squared magnitude 3 on each site.
        Then, on 1/2 the sites, we gather and sum the eight neighboring
        random vectors and add 2*m times the local vector.
            source = M_adjoint*R, on even sites
   reinitialize after niters iterations and try once more.
   parity=EVEN = do only even sites, parity=ODD = do odd sites,
   parity=EVENANDODD = do all sites
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/prefetch.h"
#define FETCH_UP 1
/* dslash_fn_tmp.c */
void cleanup_gathers(msg_tag *t1[16UL],msg_tag *t2[16UL]);
#define LOOPEND
#include "../include/loopend.h"
#ifdef CONGRAD_TMP_VECTORS
su3_vector *ttt;
su3_vector *cg_p;
su3_vector *resid;
su3_vector *t_dest;
int first_congrad = 1;
#endif

int ks_congrad(field_offset src,field_offset dest,float mass,int niter,float rsqmin,int parity,float *final_rsq_ptr)
{
  register int i;
  register site *s;
/* counter for iterations */
  int iteration;
/* Sugar's a,b,resid**2,last resid*2 */
  float a;
  float b;
/* rsq = |resid|, true_rsq = rsq from actual
				   summation of resid */
  double true_rsq;
  double rsq;
  double oldrsq;
/* pkp = cg_p.K.cg_p */
  double pkp;
/* Re<resid|ttt>, <ttt|ttt> */
  double c_tr;
  double c_tt;
  double tempsum[4UL];
/* 4*mass*mass */
  float msq_x4;
/* squared magnitude of source vector */
  double source_norm;
/* stopping residual normalized by source norm */
  double rsqstop;
/* parity we are currently doing */
  int l_parity;
/* the other parity */
  int l_otherparity;
/* tags for gathers to parity and opposite */
  msg_tag *tags1[16UL];
  msg_tag *tags2[16UL];
/* 1 if dslash_special has been called */
  int special_started;
/* Timing */
#ifdef CGTIME
  double dtimed;
  double dtimec;
/*TEST*/
  double reduce_time;
#endif
  double nflop;
/* debug */
#ifdef CGTIME
  dtimec = -dclock();
/*TEST*/
  reduce_time = 0.0;
#endif
  nflop = 1187;
  if (parity == 3) 
    nflop *= 2;
  special_started = 0;
/* if we want both parities, we will do even first. */
  switch(parity){
    case 2:
{
      l_parity = 2;
      l_otherparity = 1;
      break; 
    }
    case 1:
{
      l_parity = 1;
      l_otherparity = 2;
      break; 
    }
    case 3:
{
      l_parity = 2;
      l_otherparity = 1;
      break; 
    }
  }
  msq_x4 = ((4.0 * mass) * mass);
  iteration = 0;
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
#ifdef CONGRAD_TMP_VECTORS
/* now we can allocate temporary variables and copy then */
/* PAD may be used to avoid cache trashing */
#define PAD 0
  if (first_congrad != 0) {
    ttt = ((su3_vector *)(malloc(((sites_on_node + 0) * sizeof(su3_vector )))));
    cg_p = ((su3_vector *)(malloc(((sites_on_node + 0) * sizeof(su3_vector )))));
    resid = ((su3_vector *)(malloc(((sites_on_node + 0) * sizeof(su3_vector )))));
    t_dest = ((su3_vector *)(malloc(((sites_on_node + 0) * sizeof(su3_vector )))));
    first_congrad = 0;
  }
#endif
#ifdef CGTIME
  dtimec = -dclock();
#endif
#ifdef CONGRAD_TMP_VECTORS
/* now we copy dest to temporaries */
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    t_dest[i] =  *((su3_vector *)(((char *)s) + dest));
  }
#endif
/* initialization process */
  start:
/**node0_printf("ks_congrad4: start, parity = %d\n",parity);**/
/* ttt <-  (-1)*M_adjoint*M*dest
           resid,cg_p <- src + ttt
           rsq = |resid|^2
           source_norm = |src|^2
        */
/* clean up gathers */
  if (special_started == 1) {
    cleanup_gathers(tags1,tags2);
    special_started = 0;
  }
/**if(this_node==0)if(iteration>1)printf("CONGRAD: restart rsq = %.10e\n",rsq);**/
  rsq = (source_norm = 0.0);
#ifdef CONGRAD_TMP_VECTORS
  dslash_fn_on_temp_special(t_dest,ttt,l_otherparity,tags2,1);
  dslash_fn_on_temp_special(ttt,ttt,l_parity,tags1,1);
  cleanup_gathers(tags1,tags2);
#else
/* Why not use dslash_fn_special here ??? -CD */
#endif
/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
{
    register int loopend;
    loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
#ifdef CONGRAD_TMP_VECTORS
      if (i < (loopend - 1)) {;
      }
      scalar_mult_add_su3_vector((ttt + i),(t_dest + i),-msq_x4,(ttt + i));
/* note that we go back to the site structure for src */
      add_su3_vector(((su3_vector *)(((char *)s) + src)),(ttt + i),(resid + i));
/* remember ttt contains -M_adjoint*M*src */
      cg_p[i] = resid[i];
/* note that we go back to the site structure for src */
      source_norm += ((double )(magsq_su3vec(((su3_vector *)(((char *)s) + src)))));
      rsq += ((double )(magsq_su3vec((resid + i))));
#else
#endif
    }
  }
#ifdef CGTIME
  reduce_time -= dclock();
#endif
/* not yet summed over nodes */
  true_rsq = rsq;
  tempsum[0] = source_norm;
  tempsum[1] = rsq;
  g_vecdoublesum(tempsum,2);
  source_norm = tempsum[0];
  rsq = tempsum[1];
#ifdef CGTIME
  reduce_time += dclock();
#endif
/**if(this_node==0)printf("CONGRAD: start rsq = %.10e\n",rsq);**/
/* iteration counts number of multiplications
                           by M_adjoint*M */
  iteration++;
  total_iters++;
  rsqstop = (rsqmin * source_norm);
/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/
  if (rsq <= rsqstop) {
/* if parity==EVENANDODD, set up to do odd sites and go back */
    if (parity == 3) {
      l_parity = 1;
      l_otherparity = 2;
/* so we won't loop endlessly */
      parity = 2;
      iteration = 0;
/**node0_printf("instant goto start\n"); **/
      goto start;
    }
     *final_rsq_ptr = ((float )rsq);
/**node0_printf("instant return\n"); fflush(stdout);**/
    return iteration;
  }
/**pkp=0.0;
	if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	iteration,(double)rsq,(double)pkp);fflush(stdout);}**/
/* main loop - do until convergence or time to restart */
/*
           oldrsq <- rsq
           ttt <- (-1)*M_adjoint*M*cg_p
           pkp <- (-1)*cg_p.M_adjoint*M.cg_p
           a <- -rsq/pkp
           dest <- dest + a*cg_p
           resid <- resid + a*ttt
           rsq <- |resid|^2
           b <- rsq/oldrsq
           cg_p <- resid + b*cg_p
        */
  do {
/* not yet summed over nodes */
    oldrsq = true_rsq;
    pkp = 0.0;
/* sum of neighbors */
    if (special_started == 0) {
#ifdef CONGRAD_TMP_VECTORS
      dslash_fn_on_temp_special(cg_p,ttt,l_otherparity,tags2,1);
      dslash_fn_on_temp_special(ttt,ttt,l_parity,tags1,1);
#else
#endif
      special_started = 1;
    }
    else {
#ifdef CONGRAD_TMP_VECTORS
      dslash_fn_on_temp_special(cg_p,ttt,l_otherparity,tags2,0);
      dslash_fn_on_temp_special(ttt,ttt,l_parity,tags1,0);
#else
#endif
    }
/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
/* pkp  <- cg_p.(ttt - msq*cg_p) */
    pkp = 0.0;
    c_tr = 0.0;
    c_tt = 0.0;
{
      register int loopend;
      loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
      for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
#ifdef CONGRAD_TMP_VECTORS
        if (i < (loopend - 1)) {;
        }
        scalar_mult_add_su3_vector((ttt + i),(cg_p + i),-msq_x4,(ttt + i));
        pkp += ((double )(su3_rdot((cg_p + i),(ttt + i))));
        c_tr += ((double )(su3_rdot((ttt + i),(resid + i))));
        c_tt += ((double )(su3_rdot((ttt + i),(ttt + i))));
#else
#endif
      }
    }
#ifdef CGTIME
    reduce_time -= dclock();
#endif
/* finally sum oldrsq over nodes, also other sums */
    tempsum[0] = pkp;
    tempsum[1] = c_tr;
    tempsum[2] = c_tt;
    tempsum[3] = oldrsq;
    g_vecdoublesum(tempsum,4);
    pkp = tempsum[0];
    c_tr = tempsum[1];
    c_tt = tempsum[2];
    oldrsq = tempsum[3];
#ifdef CGTIME
    reduce_time += dclock();
#endif
    iteration++;
    total_iters++;
    a = ((float )(-rsq / pkp));
/* dest <- dest - a*cg_p */
/* resid <- resid - a*ttt */
    true_rsq = 0.0;
{
      register int loopend;
      loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
      for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
#ifdef CONGRAD_TMP_VECTORS
        if (i < (loopend - 1)) {;
        }
        scalar_mult_add_su3_vector((t_dest + i),(cg_p + i),a,(t_dest + i));
        scalar_mult_add_su3_vector((resid + i),(ttt + i),a,(resid + i));
        true_rsq += ((double )(magsq_su3vec((resid + i))));
#else
#endif
      }
    }
/**printf("XXX:  node %d\t%e\t%e\t%e\n",this_node,oldrsq,c_tr,c_tt);**/
/*TEST - should equal true_rsq */
    rsq = ((oldrsq + ((2.0 * a) * c_tr)) + ((a * a) * c_tt));
/**c_tt = true_rsq;**/
/* TEMP for test */
/**g_doublesum(&c_tt);**/
/* TEMP true value for rsq */
/**node0_printf("RSQTEST: %e\t%e\t%e\n",rsq,c_tt,rsq-c_tt);**/
/**if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}**/
    if (rsq <= rsqstop) {
#ifdef CONGRAD_TMP_VECTORS
/* copy t_dest back to site structure */
{
        register int loopend;
        loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
        for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
           *((su3_vector *)(((char *)s) + dest)) = t_dest[i];
        }
      }
#endif
/* if parity==EVENANDODD, set up to do odd sites and go back */
      if (parity == 3) {
        l_parity = 1;
        l_otherparity = 2;
/* so we won't loop endlessly */
        parity = 2;
        iteration = 0;
/**node0_printf("normal goto start\n"); **/
        goto start;
      }
       *final_rsq_ptr = ((float )rsq);
      if (special_started == 1) {
        cleanup_gathers(tags1,tags2);
        special_started = 0;
      }
/**node0_printf("normal return\n"); fflush(stdout);**/
#ifdef CGTIME
      dtimec += dclock();
      if (this_node == 0) {
        printf("CONGRAD5: time = %e iters = %d mflops = %e\n",dtimec,iteration,(((nflop * volume) * iteration) / ((1.0e6 * dtimec) * (numnodes()))));
        printf("TESTCONG: reduce_time = %e iters = %d time/iter = %e\n",reduce_time,iteration,(reduce_time / iteration));
//{ /* time stamp for NERSC performance studies */
//      time_t time_now;
//      char time_out[26];
//      time(&time_now);
//      ctime_r(&time_now,time_out);
//      printf("  Time stamp:   %s",time_out);
//}
        fflush(stdout);
      }
#endif
      return iteration;
    }
    b = (((float )rsq) / oldrsq);
/* cg_p  <- resid + b*cg_p */
#ifdef CONGRAD_TMP_VECTORS
{
      register int loopend;
      loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
      for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
        scalar_mult_add_su3_vector((resid + i),(cg_p + i),b,(cg_p + i));
      }
    }
#else
#endif
  }while ((iteration % niter) != 0);
  if (iteration < (5 * niter)) {
/**node0_printf("try again goto start\n");**/
    goto start;
  }
#ifdef CONGRAD_TMP_VECTORS
/* if we have gotten here, no convergence after several restarts: must
	copy t_dest back to site structure */
{
    register int loopend;
    loopend = ((l_parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((l_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
       *((su3_vector *)(((char *)s) + dest)) = t_dest[i];
    }
  }
#endif
/* if parity==EVENANDODD, set up to do odd sites and go back */
  if (parity == 3) {
    l_parity = 1;
    l_otherparity = 2;
/* so we won't loop endlessly */
    parity = 2;
    iteration = 0;
    goto start;
  }
   *final_rsq_ptr = rsq;
/* clean up gathers */
  if (special_started == 1) {
    cleanup_gathers(tags1,tags2);
    special_started = 0;
  }
  if (this_node == 0) 
    printf("CG not converged after %d iterations, res. = %e wanted %e\n",iteration,rsq,rsqstop);
  fflush(stdout);
  return iteration;
}
/* clear an su3_vector in the lattice */

void clear_latvec(field_offset v,int parity)
{
  register int i;
  register int j;
  register site *s;
  register su3_vector *vv;
  switch(parity){
    case 0x02:
{
      for (((i = 0) , (s = lattice)); i < even_sites_on_node; (i++ , s++)) {
        vv = ((su3_vector *)(((char *)s) + v));
        for (j = 0; j < 3; j++) {
          (vv -> c)[j].real = ((vv -> c)[j].imag = 0.0);
        }
      }
      break; 
    }
    case 0x01:
{
      for (((i = even_sites_on_node) , (s = (lattice + i))); i < sites_on_node; (i++ , s++)) {
        vv = ((su3_vector *)(((char *)s) + v));
        for (j = 0; j < 3; j++) {
          (vv -> c)[j].real = ((vv -> c)[j].imag = 0.0);
        }
      }
      break; 
    }
    case 0x03:
{
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        vv = ((su3_vector *)(((char *)s) + v));
        for (j = 0; j < 3; j++) {
          (vv -> c)[j].real = ((vv -> c)[j].imag = 0.0);
        }
      }
      break; 
    }
  }
}
/* copy an su3_vector in the lattice */

void copy_latvec(field_offset src,field_offset dest,int parity)
{
  register int i;
  register site *s;
  register su3_vector *spt;
  register su3_vector *dpt;
  switch(parity){
    case 0x02:
{
      for (((i = 0) , (s = lattice)); i < even_sites_on_node; (i++ , s++)) {
        s = (lattice + i);
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
         *dpt =  *spt;
      }
      break; 
    }
    case 0x01:
{
      for (((i = even_sites_on_node) , (s = (lattice + i))); i < sites_on_node; (i++ , s++)) {
        s = (lattice + i);
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
         *dpt =  *spt;
      }
      break; 
    }
    case 0x03:
{
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        s = (lattice + i);
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
         *dpt =  *spt;
      }
      break; 
    }
  }
}
/* scalar multiply and add an SU3 vector in the lattice */

void scalar_mult_add_latvec(field_offset src1,field_offset src2,float scalar,field_offset dest,int parity)
{
  register int i;
  register site *s;
  register su3_vector *spt1;
  register su3_vector *spt2;
  register su3_vector *dpt;
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      spt1 = ((su3_vector *)(((char *)s) + src1));
      spt2 = ((su3_vector *)(((char *)s) + src2));
      dpt = ((su3_vector *)(((char *)s) + dest));
      if (i < (loopend - 1)) {;
      }
      scalar_mult_add_su3_vector(spt1,spt2,scalar,dpt);
    }
  }
}

void scalar2_mult_add_su3_vector(su3_vector *a,float s1,su3_vector *b,float s2,su3_vector *c)
{
  register int i;
  for (i = 0; i < 3; i++) {
    (c -> c)[i].real = ((s1 * (a -> c)[i].real) + (s2 * (b -> c)[i].real));
    (c -> c)[i].imag = ((s1 * (a -> c)[i].imag) + (s2 * (b -> c)[i].imag));
  }
}
/* scalar multiply two SU3 vector and add through the lattice */

void scalar2_mult_add_latvec(field_offset src1,float scalar1,field_offset src2,float scalar2,field_offset dest,int parity)
{
  register int i;
  register site *s;
  register su3_vector *spt1;
  register su3_vector *spt2;
  register su3_vector *dpt;
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      spt1 = ((su3_vector *)(((char *)s) + src1));
      spt2 = ((su3_vector *)(((char *)s) + src2));
      dpt = ((su3_vector *)(((char *)s) + dest));
      if (i < (loopend - 1)) {;
      }
      scalar2_mult_add_su3_vector(spt1,scalar1,spt2,scalar2,dpt);
    }
  }
}
/* scalar multiply an SU3 vector in the lattice */

void scalar_mult_latvec(field_offset src,float scalar,field_offset dest,int parity)
{
  register int i;
  register site *s;
  register su3_vector *spt;
  register su3_vector *dpt;
  switch(parity){
    case 0x02:
{
      for (((i = 0) , (s = lattice)); i < even_sites_on_node; (i++ , s++)) {
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
        scalar_mult_su3_vector(spt,scalar,dpt);
      }
      break; 
    }
    case 0x01:
{
      for (((i = even_sites_on_node) , (s = (lattice + i))); i < sites_on_node; (i++ , s++)) {
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
        scalar_mult_su3_vector(spt,scalar,dpt);
      }
      break; 
    }
    case 0x03:
{
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        spt = ((su3_vector *)(((char *)s) + src));
        dpt = ((su3_vector *)(((char *)s) + dest));
        scalar_mult_su3_vector(spt,scalar,dpt);
      }
      break; 
    }
  }
}
