/******* d_congrad5_fn_tmp.c - conjugate gradient for SU3/fermions ****/
/* TEST VERSION  4/18/03, TIMING FOR GLOBAL REDUCTIONS */
/* TRY HOMEBREW GLOBAL SUM */
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

/*TEST*/
#include <mpi.h>
#define SEND_FIELD_ID      2  /* id of field sent from one node to another */
void g_doublesum_test( double *dpt );
/* for the moment, try to be efficient for numbers of processors
 * that are powers of 2 time some small number 
*/
void g_doublesum_test( double *dpt ){
    double sum, buf, *sumlist; /* sum, and array of partial sums */
    int i, mask, nodes, partner, dim, mypos, tonode, fromnode;
    MPI_Status status;

    nodes = numnodes();
    sum = *dpt;
    /* first, pairwise sums, for powers of two in num_nodes */
    for( mask=0x1; nodes & mask == 0; mask <<= 1 ){
	partner = this_node ^ mask;
	buf = sum;
	MPI_Sendrecv_replace( &buf, 1, MPI_DOUBLE, partner, SEND_FIELD_ID,
	    partner, SEND_FIELD_ID, MPI_COMM_WORLD, &status );
	sum += buf;

    }

    /* sum over remaining "directions".  Remember addition not associative */
    dim = nodes/mask;	/* number of subsets left to be summed */
    mypos = (this_node + mask - 1 )/mask; /* this processors position in list */
    sumlist = (double *)malloc( dim*sizeof(double) );
    /* pass array around cyclically, fill in my contribution each time */
    tonode = (this_node + mask) % nodes;
    fromnode = (this_node + nodes - mask ) % nodes;
    for( i=0; i<dim ; i++ ){
	sumlist[mypos] = sum;
	 MPI_Sendrecv_replace( sumlist, dim, MPI_DOUBLE, tonode, SEND_FIELD_ID,
	    fromnode, SEND_FIELD_ID, MPI_COMM_WORLD, &status );
    }
    /* add.  all processors must sum in same order */
    for( i=0,sum=0.0; i<dim; i++)sum += sumlist[i];

    free(sumlist);
    *dpt = sum;
}
/*END TEST*/


void cleanup_gathers(msg_tag *t1[16],msg_tag *t2[16]); /* dslash_fn_tmp.c */

#define LOOPEND
#include "../include/loopend.h"

#ifdef CONGRAD_TMP_VECTORS
su3_vector *ttt,*cg_p;
su3_vector *resid;
su3_vector *t_dest;
int first_congrad = 1;
#endif

int ks_congrad( field_offset src, field_offset dest, float mass,
    int niter, float rsqmin, int parity, float *final_rsq_ptr ){
    register int i;
    register site *s;
    int iteration;	/* counter for iterations */
    float a,b;			/* Sugar's a,b,resid**2,last resid*2 */
    double rsq,oldrsq,pkp;		/* pkp = cg_p.K.cg_p */
    float msq_x4;	/* 4*mass*mass */
    double source_norm;	/* squared magnitude of source vector */
    double rsqstop;	/* stopping residual normalized by source norm */
    int l_parity;	/* parity we are currently doing */
    int l_otherparity;	/* the other parity */
    msg_tag * tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
    int special_started;	/* 1 if dslash_special has been called */

/* Timing */

#ifdef CGTIME
double dtimed,dtimec;
double reduce_time,reduce_newtime; /*TEST*/
#endif
double nflop;

/* debug */
#ifdef CGTIME
 dtimec = -dclock(); 
 reduce_time = 0.0; /*TEST*/
 reduce_newtime = 0.0; /*TEST*/
#endif

nflop = 1187;
if(parity==EVENANDODD)nflop *=2;
	
	special_started=0;
	/* if we want both parities, we will do even first. */
	switch(parity){
	    case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	    case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	    case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
	}
	msq_x4 = 4.0*mass*mass;
        iteration = 0;

        if (!valid_longlinks) load_longlinks();
        if (!valid_fatlinks) load_fatlinks();
#ifdef CONGRAD_TMP_VECTORS
	/* now we can allocate temporary variables and copy then */
	/* PAD may be used to avoid cache trashing */
#define PAD 0

 	if(first_congrad) {
	  ttt = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  cg_p = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  resid = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  t_dest = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  first_congrad = 0;
 	}
#endif

#ifdef CGTIME
 dtimec = -dclock(); 
#endif

#ifdef CONGRAD_TMP_VECTORS
	/* now we copy dest to temporaries */
  FORALLSITES(i,s) {
    t_dest[i] = *(su3_vector *)F_PT(s,dest);
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
	if(special_started==1) {	/* clean up gathers */
	    cleanup_gathers(tags1,tags2);
	    special_started=0;
	}
/**if(this_node==0)if(iteration>1)printf("CONGRAD: restart rsq = %.10e\n",rsq);**/
        rsq = source_norm = 0.0;
#ifdef CONGRAD_TMP_VECTORS
	dslash_fn_on_temp_special(t_dest, ttt,l_otherparity,tags2,1);
	dslash_fn_on_temp_special(ttt,ttt,l_parity,tags1,1);
	cleanup_gathers(tags1,tags2);
#else
	/* Why not use dslash_fn_special here ??? -CD */
	dslash_fn( dest, F_OFFSET(ttt), l_otherparity);
	dslash_fn(F_OFFSET(ttt),F_OFFSET(ttt),l_parity);
#endif
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITY(i,s,l_parity){
#ifdef CONGRAD_TMP_VECTORS
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &ttt[i+FETCH_UP], 
			   &t_dest[i+FETCH_UP],
			   (su3_vector *)F_PT(s+FETCH_UP,src),
			   &resid[i+FETCH_UP]);
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &t_dest[i],
				      -msq_x4, &ttt[i] );
	    /* note that we go back to the site structure for src */
	  add_su3_vector( (su3_vector *)F_PT(s,src),
			  &ttt[i], &resid[i] );
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  /* note that we go back to the site structure for src */
	  source_norm += (double)magsq_su3vec( (su3_vector *)F_PT(s,src) );
	  rsq += (double)magsq_su3vec( &resid[i] );
#else
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &((s+FETCH_UP)->ttt), 
			   (su3_vector *)F_PT(s+FETCH_UP,dest),
			   (su3_vector *)F_PT(s+FETCH_UP,src),
			   &((s+FETCH_UP)->resid));
	  }
	  scalar_mult_add_su3_vector( &(s->ttt), (su3_vector *)F_PT(s,dest),
				      -msq_x4, &(s->ttt) );
	  add_su3_vector( (su3_vector *)F_PT(s,src),
			  &(s->ttt), &(s->resid) );
	  s->cg_p = s->resid;
	  source_norm += (double) magsq_su3vec( (su3_vector *)F_PT(s,src) );
	  rsq += (double) magsq_su3vec( &(s->resid) );
#endif
	} END_LOOP
reduce_time -= dclock();
	g_doublesum( &source_norm );
        g_doublesum( &rsq );
reduce_time += dclock();
	/**if(this_node==0)printf("CONGRAD: start rsq = %.10e\n",rsq);**/
        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/
        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**node0_printf("instant goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(float)rsq;
	    /**node0_printf("instant return\n"); fflush(stdout);**/
             return (iteration);
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
    do{
        oldrsq = rsq;
        pkp = 0.0;
	/* sum of neighbors */

	if(special_started==0){
#ifdef CONGRAD_TMP_VECTORS
	    dslash_fn_on_temp_special( cg_p, ttt, l_otherparity, tags2, 1 );
	    dslash_fn_on_temp_special( ttt, ttt, l_parity, tags1, 1);
#else
	    dslash_fn_special( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity,
		tags2, 1 );
	    dslash_fn_special( F_OFFSET(ttt), F_OFFSET(ttt), l_parity,
		tags1, 1);
#endif
	    special_started=1;
	}
	else {
#ifdef CONGRAD_TMP_VECTORS
	    dslash_fn_on_temp_special( cg_p, ttt, l_otherparity, tags2, 0 );
	    dslash_fn_on_temp_special( ttt, ttt, l_parity, tags1, 0);
#else
	    dslash_fn_special( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity,
		tags2, 0 );
	    dslash_fn_special( F_OFFSET(ttt), F_OFFSET(ttt), l_parity,
		tags1, 0 );
#endif
	}

	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
#ifdef CONGRAD_TMP_VECTORS
	  if( i < loopend-FETCH_UP ){
	    prefetch_VV( &ttt[i+FETCH_UP], &cg_p[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &cg_p[i], -msq_x4,
				      &ttt[i] );
	  pkp += (double)su3_rdot( &cg_p[i], &ttt[i] );
#else
	  if( i < loopend-FETCH_UP ){
	    prefetch_VV( &((s+FETCH_UP)->ttt), &((s+FETCH_UP)->cg_p) );
	  }
	  scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4,
				      &(s->ttt) );
	  pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
#endif
	} END_LOOP
/*TEST*/
{double pkptmp;
pkptmp = pkp;
reduce_time -= dclock();
	g_doublesum( &pkp );
reduce_time += dclock();
reduce_newtime -= dclock();
g_doublesum_test( &pkptmp );
reduce_newtime += dclock();
if(this_node==  0)printf("SUMTEST_000:  %.14le\t%.14le\t%e\n",pkp,pkptmp,pkp-pkptmp);
if(this_node== 77)printf("SUMTEST_077:  %.14le\t%.14le\t%e\n",pkp,pkptmp,pkp-pkptmp);
/*END TEST*/
}
	iteration++;
	total_iters++;

	a = (float) (-rsq/pkp);

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
	rsq=0.0;
	FORSOMEPARITY(i,s,l_parity){
#ifdef CONGRAD_TMP_VECTORS
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &t_dest[i+FETCH_UP], 
			   &cg_p[i+FETCH_UP], 
			   &resid[i+FETCH_UP], 
			   &ttt[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
	  scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
	  rsq += (double)magsq_su3vec( &resid[i] );
#else
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( (su3_vector *)F_PT((s+FETCH_UP),dest), 
			   &((s+FETCH_UP)->cg_p),
			   &((s+FETCH_UP)->resid), 
			   &((s+FETCH_UP)->ttt) );
	  }
	  scalar_mult_add_su3_vector( (su3_vector *)F_PT(s,dest),
				      &(s->cg_p), a, (su3_vector *)F_PT(s,dest) );
	  scalar_mult_add_su3_vector( &(s->resid), &(s->ttt), a, &(s->resid));
	  rsq += (double)magsq_su3vec( &(s->resid) );
#endif
	} END_LOOP
{/*TEST*/
double rsqtmp;
rsqtmp=rsq;
reduce_time -= dclock();
	g_doublesum(&rsq);
reduce_time += dclock();
reduce_newtime -= dclock();
	g_doublesum_test(&rsqtmp);
reduce_newtime += dclock();
}/*ENDTEST*/
	/**if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}**/
	
        if( rsq <= rsqstop ){
#ifdef CONGRAD_TMP_VECTORS
	  /* copy t_dest back to site structure */
          FORSOMEPARITY(i,s,l_parity){
                  *(su3_vector *)F_PT(s,dest) = t_dest[i];
          } END_LOOP
#endif
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**node0_printf("normal goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(float)rsq;
	    if(special_started==1) {
	      cleanup_gathers(tags1,tags2);
	      special_started = 0;
	    }
	    
	    /**node0_printf("normal return\n"); fflush(stdout);**/
#ifdef CGTIME
 dtimec += dclock();
if(this_node==0){printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
dtimec,iteration,(double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
printf("TESTCONG: reduce_time = %e iters = %d time/iter = %e\n",
reduce_time,iteration,reduce_time/iteration );
printf("TESTCONG: reduce_newtime = %e iters = %d time/iter = %e\n",
reduce_newtime,iteration,reduce_newtime/iteration );
fflush(stdout);}
#endif
             return (iteration);
        }

	b = (float)rsq/oldrsq;
	/* cg_p  <- resid + b*cg_p */
#ifdef CONGRAD_TMP_VECTORS
        FORSOMEPARITY(i,s,l_parity){
           scalar_mult_add_su3_vector( &resid[i],
                                      &cg_p[i] , b , &cg_p[i]);
        } END_LOOP
#else
	scalar_mult_add_latvec( F_OFFSET(resid), F_OFFSET(cg_p),
	    b, F_OFFSET(cg_p), l_parity);
#endif

    } while( iteration%niter != 0);

    if( iteration < 5*niter ){
	/**node0_printf("try again goto start\n");**/
	 goto start;
    }
#ifdef CONGRAD_TMP_VECTORS
    /* if we have gotten here, no convergence after several restarts: must
	copy t_dest back to site structure */
          FORSOMEPARITY(i,s,l_parity){
                  *(su3_vector *)F_PT(s,dest) = t_dest[i];
          } END_LOOP
#endif

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=rsq;
    if(special_started==1){	/* clean up gathers */
      cleanup_gathers(tags1,tags2);
      special_started = 0;
    }
    node0_printf(
        "CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    return(iteration);
}

/* clear an su3_vector in the lattice */
void clear_latvec(field_offset v,int parity){
register int i,j;
register site *s;
register su3_vector *vv;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case ODD: FORODDSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
    } 
}

/* copy an su3_vector in the lattice */
void copy_latvec(field_offset src,field_offset dest,int parity){
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case ODD: FORODDSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
    } 
}

/* scalar multiply and add an SU3 vector in the lattice */
void scalar_mult_add_latvec(field_offset src1,field_offset src2,
			    float scalar,field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
               spt1 = (su3_vector *)F_PT(s,src1);
                spt2 = (su3_vector *)F_PT(s,src2);
                dpt = (su3_vector *)F_PT(s,dest);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV( (su3_vector *)F_PT((s+FETCH_UP),src1),
				(su3_vector *)F_PT((s+FETCH_UP),src2),
				(su3_vector *)F_PT((s+FETCH_UP),dest) );
		}
                scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt);
       } END_LOOP
}


void scalar2_mult_add_su3_vector(su3_vector *a, float s1, su3_vector *b, 
				 float s2, su3_vector *c){
register int i;
    for(i=0;i<3;i++){
        c->c[i].real = s1*a->c[i].real + s2*b->c[i].real;
        c->c[i].imag = s1*a->c[i].imag + s2*b->c[i].imag;
    }
}

/* scalar multiply two SU3 vector and add through the lattice */
void scalar2_mult_add_latvec(field_offset src1,float scalar1,
			     field_offset src2,float scalar2,
			     field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt  = (su3_vector *)F_PT(s,dest);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV((su3_vector *)F_PT((s+FETCH_UP),src1),
			       (su3_vector *)F_PT((s+FETCH_UP),src2),
			       (su3_vector *)F_PT((s+FETCH_UP),dest) );
		}
		scalar2_mult_add_su3_vector( spt1, scalar1, spt2, scalar2, dpt);
       } END_LOOP
}

/* scalar multiply an SU3 vector in the lattice */
void scalar_mult_latvec(field_offset src,float scalar,
			field_offset dest,int parity)
{
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case ODD: FORODDSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
    } 
}

