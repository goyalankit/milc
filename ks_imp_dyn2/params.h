#ifndef _PARAMS_H
#define _PARAMS_H

#include "../include/macros.h"  /* For MAXFILENAME */
#include "defines.h"

/* structure for passing simulation parameters to each node */
typedef struct {
	int stopflag;   /* 1 if it is time to stop */
    /* INITIALIZATION PARAMETERS */
	int nx,ny,nz,nt;  /* lattice dimensions */
	int iseed;	/* for random numbers */
	int nflavors1;	/* the number of flavors of first type */
	int nflavors2;	/* the number of flavors of second type */
    /*  REPEATING BLOCK */
	int warms;	/* the number of warmup trajectories */
	int trajecs;	/* the number of real trajectories */
	int steps;	/* number of steps for updating */
	int propinterval;     /* number of trajectories between measurements */
	float beta,mass1,mass2; /* gauge coupling, quark masses */
	float u0; /* tadpole parameter */
	int niter; 	/* maximum number of c.g. iterations */
	float rsqmin,rsqprop;  /* for deciding on convergence */
	float epsilon;	/* time step */
        char spectrum_request[MAX_SPECTRUM_REQUEST];   /* request list for spectral measurements */
        int source_start, source_inc, n_sources; /* source time and increment */
        int spectrum_multimom_nmasses; 
        float spectrum_multimom_low_mass;
        float spectrum_multimom_mass_step;
        int fpi_nmasses;
        float fpi_mass[MAX_FPI_NMASSES];
	int startflag;  /* what to do for beginning lattice */
	int saveflag;   /* what to do with lattice at end */
	char startfile[MAXFILENAME],savefile[MAXFILENAME];
}  params;

#endif /* _PARAMS_H */
