void indigo__create_map();
/************************* control.c *******************************/
/* MIMD version 6 */
/* Main procedure for SU3 with dynamical staggered fermions        */
/* general quark action, general gauge action */
/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of 
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */
#define CONTROL
#include "ks_imp_includes.h"	/* definitions files and prototypes */
/* Input gauge field header */
gauge_header start_lat_hdr;
char *szInputFile = "properties.in";
FILE *fpInput;

int main(int argc,char **argv)
{
  int i;
  int meascount;
  int traj_done;
  int prompt;
  float ssplaq;
  float stplaq;
  int s_iters;
  int avs_iters;
  int avspect_iters;
  int avbcorr_iters;
  double dtime;
  double dclock();
  indigo__init_(1,1);
  indigo__create_map();
  initialize_machine(argc,argv);
  g_sync();
/* set up */
  fpInput = fopen(szInputFile,"r");
  if (!(fpInput != 0)) {
    fprintf(stderr,"Error opening file...\n");
    return 1;
  }
  prompt = setup();
/* loop over input sets */
  while(readin(prompt) == 0){
/* perform warmup trajectories */
    dtime = -dclock();
    for (traj_done = 0; traj_done < warms; traj_done++) {
      update();
    }
    if (this_node == 0) 
      printf("WARMUPS COMPLETED\n");
    fflush(stdout);
/* perform measuring trajectories, reunitarizing and measuring 	*/
/* number of measurements 		*/
    meascount = 0;
    avspect_iters = (avs_iters = (avbcorr_iters = 0));
    for (traj_done = 0; traj_done < trajecs; traj_done++) {
/* do the trajectories */
      s_iters = update();
/* measure every "propinterval" trajectories */
      if ((traj_done % propinterval) == (propinterval - 1)) {
/* call gauge_variable fermion_variable measuring routines */
/* results are printed in output file */
        rephase(0);
        g_measure();
        rephase(1);
        f_meas_imp(((field_offset )(((char *)(&lattice[0].phi1)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),mass1);
        f_meas_imp(((field_offset )(((char *)(&lattice[0].phi2)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))),mass2);
#ifdef SPECTRUM 
/* Fix TUP Coulomb gauge - gauge links only*/
#ifdef HYBRIDS
#endif
#ifdef FN
#endif
#endif /* SPECTRUM */
        avs_iters += s_iters;
        ++meascount;
        fflush(stdout);
      }
/* end loop over trajectories */
    }
    if (this_node == 0) 
      printf("RUNNING COMPLETED\n");
    fflush(stdout);
    if (meascount > 0) {
      if (this_node == 0) 
        printf("average cg iters for step= %e\n",(((double )avs_iters) / meascount));
#ifdef SPECTRUM
#endif
    }
    dtime += dclock();
    if (this_node == 0) {
      printf("total_iters = %d\n",total_iters);
      printf("NERSC_TIME %8.3f secs\n",dtime);
    }
    fflush(stdout);
/* save lattice if requested */
    if (saveflag != 20) {
      rephase(0);
      save_lattice(saveflag,savefile);
      rephase(1);
    }
/*TEMP*/
#ifdef NONSENSE
#endif
  }
  fclose(fpInput);
  return 0;
}
