/******* ks_multicg.c - multi-mass CG for SU3/fermions ****/
/* MIMD version 6 */

/* Multi-mass CG inverter for staggered fermions */

/* Based on B. Jegerlehner, hep-lat/9612014.
   See also A. Frommer, S. G\"usken, T. Lippert, B. N\"ockel,"
   K. Schilling, Int. J. Mod. Phys. C6 (1995) 627. */

/* This version is based on d_congrad5_fn.c and d_congrad5_eo.c */

/* For "fat link actions", ie when FN is defined, this version
   assumes connection to nearest neighbor points is stored in fatlink.
   For actions with a Naik term, it assumes the connection to third
   nearest neighbors is in longlink. */


#include "generic_ks_includes.h"	/* definitions files and prototypes */

#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif

#include "../include/loopend.h"

int ks_multicg(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    float *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    float rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    float *final_rsq_ptr	/* final residue squared */
    )
{
    /* Site su3_vector's resid, cg_p and ttt are used as temporaies */
    register int i;
    register site *s;
    int iteration;	/* counter for iterations */
    double c1, c2, rsq, oldrsq, pkp;		/* pkp = cg_p.K.cg_p */
    double source_norm;	/* squared magnitude of source vector */
    double rsqstop;	/* stopping residual normalized by source norm */
    int l_parity;	/* parity we are currently doing */
    int l_otherparity;	/* the other parity */
    msg_tag *tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
    int special_started;	/* 1 if dslash_special has been called */
    int j, j_low;
    float *shifts, mass_low, msq_xm4;
    double *zeta_i, *zeta_im1, *zeta_ip1;
    double *beta_i, *beta_im1, *alpha;
    su3_vector **pm;	/* vectors not involved in gathers */

/* Timing */

#ifdef CGTIME
    double dtimed,dtimec;
#endif
    double nflop;

/* debug */
#ifdef CGTIME
    dtimec = -dclock(); 
#endif

    nflop = 1187;
    if(parity==EVENANDODD)nflop *=2;
	
    special_started = 0;
    /* if we want both parities, we will do even first. */
    switch(parity){
	case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
    }

    shifts = (float *)malloc(num_masses*sizeof(float));
    zeta_i = (double *)malloc(num_masses*sizeof(double));
    zeta_im1 = (double *)malloc(num_masses*sizeof(double));
    zeta_ip1 = (double *)malloc(num_masses*sizeof(double));
    beta_i = (double *)malloc(num_masses*sizeof(double));
    beta_im1 = (double *)malloc(num_masses*sizeof(double));
    alpha = (double *)malloc(num_masses*sizeof(double));

    pm = (su3_vector **)malloc(num_masses*sizeof(su3_vector *));
    mass_low = 1.0e+20;
    j_low = -1;
    for(j=0;j<num_masses;j++){
	shifts[j] = 4.0*masses[j]*masses[j];
	if (masses[j] < mass_low){
	    mass_low = masses[j];
	    j_low = j;
	}
    }
    for(j=0;j<num_masses;j++) if(j!=j_low){
	pm[j] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
	shifts[j] -= shifts[j_low];
    }
    msq_xm4 = -shifts[j_low];


    iteration = 0;

#ifdef FN
    if (!valid_longlinks) load_longlinks();
    if (!valid_fatlinks) load_fatlinks();
#endif

#ifdef CGTIME
    dtimec = -dclock(); 
#endif

    /* initialization process */
    start:
#ifdef FN
	if(special_started==1) {        /* clean up gathers */
	    cleanup_gathers(tags1, tags2);
	    special_started = 0;
	}
#endif
	source_norm = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    source_norm += (double) magsq_su3vec( (su3_vector *)F_PT(s,src) );
	    su3vec_copy((su3_vector *)F_PT(s,src), &(s->resid));
	    su3vec_copy(&(s->resid), &(s->cg_p));
	    clearvec(&(psim[j_low][i]));
	    for(j=0;j<num_masses;j++) if(j!=j_low){
		clearvec(&(psim[j][i]));
		su3vec_copy(&(s->resid), &(pm[j][i]));
	    }
	} END_LOOP
	g_doublesum( &source_norm );
	rsq = source_norm;

        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/

	for(j=0;j<num_masses;j++){
	    zeta_im1[j] = zeta_i[j] = 1.0;
	    beta_im1[j] = -1.0;
	    alpha[j] = 0.0;
	}

    do{
	oldrsq = rsq;
	/* sum of neighbors */

#ifdef FN
	if(special_started==0){
	    dslash_fn_special( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity,
		tags2, 1 );
	    dslash_fn_special( F_OFFSET(ttt), F_OFFSET(ttt), l_parity,
		tags1, 1);
	    special_started = 1;
	}
	else {
	    dslash_fn_special( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity,
		tags2, 0 );
	    dslash_fn_special( F_OFFSET(ttt), F_OFFSET(ttt), l_parity,
		tags1, 0 );
	}
#else
	dslash( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity);
	dslash( F_OFFSET(ttt), F_OFFSET(ttt), l_parity);
#endif

	/* finish computation of (-1)*M_adjoint*m*p and (-1)*p*M_adjoint*M*p */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p . ttt */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), msq_xm4,
		&(s->ttt) );
	    pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
	} END_LOOP
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	beta_i[j_low] = -rsq / pkp;

	zeta_ip1[j_low] = 1.0;
	for(j=0;j<num_masses;j++) if(j!=j_low){
	    zeta_ip1[j] = zeta_i[j] * zeta_im1[j] * beta_im1[j_low];
	    c1 = beta_i[j_low] * alpha[j_low] * (zeta_im1[j]-zeta_i[j]);
	    c2 = zeta_im1[j] * beta_im1[j_low] * (1.0+shifts[j]*beta_i[j_low]);
	    zeta_ip1[j] /= c1 + c2;
	    beta_i[j] = beta_i[j_low] * zeta_ip1[j] / zeta_i[j];
	}

	/* dest <- dest + beta*cg_p */
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(psim[j_low][i]),
		&(s->cg_p), (float)beta_i[j_low], &(psim[j_low][i]));
	    for(j=0;j<num_masses;j++) if(j!=j_low){
		scalar_mult_add_su3_vector( &(psim[j][i]),
		    &(pm[j][i]), (float)beta_i[j], &(psim[j][i]));
	    }
	} END_LOOP

	/* resid <- resid + beta*ttt */
	rsq = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(s->resid), &(s->ttt),
		(float)beta_i[j_low], &(s->resid));
	    rsq += (double)magsq_su3vec( &(s->resid) );
	} END_LOOP
	g_doublesum(&rsq);

	if( rsq <= rsqstop ){
	    /* if parity==EVENANDODD, set up to do odd sites and go back */
	    if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		goto start;
	    }
	    *final_rsq_ptr = (float)rsq;

#ifdef FN
	    if(special_started==1) {
		cleanup_gathers(tags1,tags2);
		special_started = 0;
	    }
#endif

	    /* Free stuff */
	    for(j=0;j<num_masses;j++) if(j!=j_low) free(pm[j]);
	    free(pm);

	    free(zeta_i);
	    free(zeta_ip1);
	    free(zeta_im1);
	    free(beta_i);
	    free(beta_im1);
	    free(alpha);
	    free(shifts);

#ifdef CGTIME
	    dtimec += dclock();
	    if(this_node==0){printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
dtimec,iteration,(double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
		fflush(stdout);}
#endif
             return (iteration);
        }

	alpha[j_low] = rsq / oldrsq;

	for(j=0;j<num_masses;j++) if(j!=j_low){
	    alpha[j] = alpha[j_low] * zeta_ip1[j] * beta_i[j] /
		       (zeta_i[j] * beta_i[j_low]);
	}

	/* cg_p  <- resid + alpha*cg_p */
	FORSOMEPARITY(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(s->resid), &(s->cg_p),
		(float)alpha[j_low], &(s->cg_p));
	    for(j=0;j<num_masses;j++) if(j!=j_low){
		scalar_mult_su3_vector( &(s->resid),
		    (float)zeta_ip1[j], &(s->ttt));
		scalar_mult_add_su3_vector( &(s->ttt), &(pm[j][i]),
		    (float)alpha[j], &(pm[j][i]));
	    }
	} END_LOOP

	/* scroll the scalars */
	for(j=0;j<num_masses;j++){
	    beta_im1[j] = beta_i[j];
	    zeta_im1[j] = zeta_i[j];
	    zeta_i[j] = zeta_ip1[j];
	}

    } while( iteration < niter );

    node0_printf(
	"CG not converged after %d iterations, res. = %e wanted %e\n",
	iteration, rsq, rsqstop);
    fflush(stdout);

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=rsq;

#ifdef FN
    if(special_started==1){	/* clean up gathers */
	cleanup_gathers(tags1, tags2);
	special_started = 0;
    }
#endif

    /* Free stuff */
    for(j=0;j<num_masses;j++) if(j!=j_low) free(pm[j]);
    free(pm);

    free(zeta_i);
    free(zeta_ip1);
    free(zeta_im1);
    free(beta_i);
    free(beta_im1);
    free(alpha);
    free(shifts);

    return(iteration);
}

