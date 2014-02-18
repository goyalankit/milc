/********** update.c ****************************************************/
/* MIMD version 6 */
/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed
 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update()
{
  int step;
  int iters = 0;
  float final_rsq;
#ifdef HMC_ALGORITHM
#endif
/* refresh the momenta */
  ranmom();
/* do "steps" microcanonical steps"  */
  for (step = 1; step <= steps; step++) {
#ifdef PHI_ALGORITHM
/* generate a pseudofermion configuration only at start*/
/* also clear xxx, since zero is our best guess for the solution
	   with a new random phi field. */
#ifdef HMC_ALGORITHM
/* find action */
/* do conjugate gradient to get (Madj M)inverse * phi */
/* do conjugate gradient to get (Madj M)inverse * phi */
/* copy link field to old_link */
#endif
/* update U's to middle of interval */
#else /* "R" algorithm */
/* first update the U's to special time interval */
/* and generate a pseudofermion configuration */
/* probably makes most sense if nflavors1 >= nflavors2 */
    update_u((epsilon * (0.5 - (nflavors1 / 8.0))));
    clear_latvec(((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),3);
    grsource_imp(((field_offset )(((char *)(&lattice[0].phi1)) - ((char *)(lattice + 0)))),mass1,2);
    update_u((epsilon * ((nflavors1 - nflavors2) / 8.0)));
    clear_latvec(((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))),3);
    grsource_imp(((field_offset )(((char *)(&lattice[0].phi2)) - ((char *)(lattice + 0)))),mass2,2);
/* update U's to middle of interval */
    update_u(((epsilon * nflavors2) / 8.0));
#endif
/* do conjugate gradient to get (Madj M)inverse * phi */
    iters += ks_congrad(((field_offset )(((char *)(&lattice[0].phi1)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),mass1,niter,rsqmin,2,&final_rsq);
    iters += ks_congrad(((field_offset )(((char *)(&lattice[0].phi2)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))),mass2,niter,rsqmin,2,&final_rsq);
/* now update H by full time interval */
    update_h(epsilon);
/* update U's by half time step to get to even time */
    update_u((epsilon * 0.5));
/* reunitarize the gauge field */
    rephase(0);
    reunitarize();
    rephase(1);
/* end loop over microcanonical steps */
  }
#ifdef HMC_ALGORITHM
/* find action */
/* do conjugate gradient to get (Madj M)inverse * phi */
/* decide whether to accept, if not, copy old link field back */
/* careful - must generate only one random number for whole lattice */
#ifdef FN
#endif
#endif
  if (steps > 0) 
    return iters / steps;
  else 
    return -99;
}
