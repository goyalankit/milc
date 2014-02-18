/************************ setup.c ****************************/
/* MIMD version 6 */
/*			    -*- Mode: C -*-
// File: setup.c
// Created: Fri Aug  4 1995
// Authors: J. Hetrick & K. Rummukainen
// Modified for general improved action 5/24/97  DT
//
// Description: Setup routines for improved fermion lattices
//              Includes lattice structures for Naik imroved 
//              staggered Dirac operator
//         Ref: S. Naik, Nucl. Phys. B316 (1989) 238
//              Includes a parameter prompt for Lepage-Mackenzie 
//              tadpole improvement
//         Ref: Phys. Rev. D48 (1993) 2250
//
*/
/* MIMD version 6 */
#define IF_OK if(status==0)
#include "ks_imp_includes.h"	/* definitions files and prototypes */
extern gauge_header start_lat_hdr;
gauge_file *gf;
gauge_file *r_parallel_i(char *);
void r_parallel(gauge_file *,field_offset );
void r_parallel_f(gauge_file *);
gauge_file *r_binary_i(char *);
void r_binary(gauge_file *);
void r_binary_f(gauge_file *);
/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int setup()
{
  int initial_set();
  void make_3n_gathers();
  int prompt;
/* print banner, get volume, nflavors1,nflavors2, seed */
  prompt = initial_set();
/* initialize the node random number generator */
  initialize_prn(&node_prn,iseed,(volume + mynode()));
/* Initialize the layout functions, which decide where sites live */
  setup_layout();
/* allocate space for lattice, set up coordinate fields */
  make_lattice();
  if (this_node == 0) 
    printf("Made lattice\n");
  fflush(stdout);
/* set up neighbor pointers and comlink structures
	   code for this routine is in com_machine.c  */
  make_nn_gathers();
  if (this_node == 0) 
    printf("Made nn gathers\n");
  fflush(stdout);
/* set up 3rd nearest neighbor pointers and comlink structures
	   code for this routine is below  */
  make_3n_gathers();
  if (this_node == 0) 
    printf("Made 3nn gathers\n");
  fflush(stdout);
/* set up K-S phase vectors, boundary conditions */
  phaseset();
  if (this_node == 0) 
    printf("Finished setup\n");
  fflush(stdout);
  return prompt;
}
/* SETUP ROUTINES */

int initial_set()
{
  int prompt;
  int status;
/* On node zero, read lattice size, seed, nflavors1, nflavors2
	 and send to others */
  if (mynode() == 0) {
/* print banner */
    printf("SU3 with improved KS action\n");
    printf("Microcanonical simulation with refreshing\n");
    printf("MIMD version 6\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
#endif
#ifdef PHI_ALGORITHM
#else
    printf("R algorithm\n");
#endif
#ifdef SPECTRUM
#endif
    status = get_prompt(&prompt);
    if (status == 0) 
      status += get_i(prompt,"nflavors1",&par_buf.nflavors1);
    if (status == 0) 
      status += get_i(prompt,"nflavors2",&par_buf.nflavors2);
#ifdef PHI_ALGORITHM
#endif
    if (status == 0) 
      status += get_i(prompt,"nx",&par_buf.nx);
    if (status == 0) 
      status += get_i(prompt,"ny",&par_buf.ny);
    if (status == 0) 
      status += get_i(prompt,"nz",&par_buf.nz);
    if (status == 0) 
      status += get_i(prompt,"nt",&par_buf.nt);
    if (status == 0) 
      status += get_i(prompt,"iseed",&par_buf.iseed);
    if (status > 0) 
      par_buf.stopflag = 1;
    else 
      par_buf.stopflag = 0;
/* end if(mynode()==0) */
  }
/* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes(((char *)(&par_buf)),(sizeof(par_buf)));
  if (par_buf.stopflag != 0) 
    normal_exit(0);
  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;
  nflavors1 = par_buf.nflavors1;
  nflavors2 = par_buf.nflavors2;
  this_node = mynode();
  number_of_nodes = numnodes();
  volume = (((nx * ny) * nz) * nt);
  total_iters = 0;
  return prompt;
}
/* read in parameters and coupling constants	*/

int readin(int prompt)
{
/* read in parameters for su3 monte carlo	*/
/* argument "prompt" is 1 if prompts are to be given for input	*/
  int status;
  float x;
  int i;
  char request_buf[512UL];
/* On node zero, read parameters and send to all other nodes */
  if (this_node == 0) {
    printf("\n\n");
    status = 0;
/* warms, trajecs */
    if (status == 0) 
      status += get_i(prompt,"warms",&par_buf.warms);
    if (status == 0) 
      status += get_i(prompt,"trajecs",&par_buf.trajecs);
/* trajectories between propagator measurements */
    if (status == 0) 
      status += get_i(prompt,"traj_between_meas",&par_buf.propinterval);
/* get couplings and broadcast to nodes	*/
/* beta, mass1, mass2 */
    if (status == 0) 
      status += get_f(prompt,"beta",&par_buf.beta);
    if (status == 0) 
      status += get_f(prompt,"mass1",&par_buf.mass1);
    if (status == 0) 
      status += get_f(prompt,"mass2",&par_buf.mass2);
    if (status == 0) 
      status += get_f(prompt,"u0",&par_buf.u0);
/* microcanonical time step */
    if (status == 0) 
      status += get_f(prompt,"microcanonical_time_step",&par_buf.epsilon);
/*microcanonical steps per trajectory */
    if (status == 0) 
      status += get_i(prompt,"steps_per_trajectory",&par_buf.steps);
/* maximum no. of conjugate gradient iterations */
    if (status == 0) 
      status += get_i(prompt,"max_cg_iterations",&par_buf.niter);
/* error per site for conjugate gradient */
    if (status == 0) 
      status += get_f(prompt,"error_per_site",&x);
/* rsqmin is r**2 in conjugate gradient */
    if (status == 0) 
      par_buf.rsqmin = (x * x);
/* New conjugate gradient normalizes rsqmin by norm of source */
/* error for propagator conjugate gradient */
    if (status == 0) 
      status += get_f(prompt,"error_for_propagator",&x);
    if (status == 0) 
      par_buf.rsqprop = (x * x);
#ifdef SPECTRUM
/* request list for spectral measurments */
/* prepend and append a comma for ease in parsing */
/* source time slice and increment */
/* Additional parameters for spectrum_multimom */
/* Additional parameters for fpi */
#endif /*SPECTRUM*/
/* find out what kind of starting lattice to use */
    if (status == 0) 
      status += ask_starting_lattice(prompt,&par_buf.startflag,par_buf.startfile);
/* find out what to do with lattice at end */
    if (status == 0) 
      status += ask_ending_lattice(prompt,&par_buf.saveflag,par_buf.savefile);
    if (status > 0) 
      par_buf.stopflag = 1;
    else 
      par_buf.stopflag = 0;
/* end if(this_node==0) */
  }
/* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes(((char *)(&par_buf)),(sizeof(par_buf)));
  if (par_buf.stopflag != 0) 
    normal_exit(0);
  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  propinterval = par_buf.propinterval;
  niter = par_buf.niter;
  rsqmin = par_buf.rsqmin;
  rsqprop = par_buf.rsqprop;
  epsilon = par_buf.epsilon;
  beta = par_buf.beta;
  mass1 = par_buf.mass1;
  mass2 = par_buf.mass2;
  u0 = par_buf.u0;
#ifdef SPECTRUM
#endif /*SPECTRUM*/
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
/* Do whatever is needed to get lattice */
  if (startflag == 10) {
    rephase(0);
  }
  startlat_p = reload_lattice(startflag,startfile);
/* if a lattice was read in, put in KS phases and AP boundary condition */
  valid_fatlinks = (valid_longlinks = 0);
  phases_in = 0;
  rephase(1);
/* make table of coefficients and permutations of loops in gauge action */
  make_loop_table();
/* make table of coefficients and permutations of paths in quark action */
  make_path_table();
  return 0;
}
/* Set up comlink structures for 3rd nearest gather pattern; 
   make_lattice() and  make_nn_gathers() must be called first, 
   preferably just before calling make_3n_gathers().
 */

void make_3n_gathers()
{
  int i;
  void *tt[16UL];
  void third_neighbor(int ,int ,int ,int ,int *,int ,int *,int *,int *,int *);
  for (i = 0; i <= 3; i++) {
    make_gather(third_neighbor,&i,1,0,1);
  }
/* Sort into the order we want for nearest neighbor gathers,
       so you can use X3UP, X3DOWN, etc. as argument in calling them. */
  sort_eight_gathers(8);
}
/* this routine uses only fundamental directions (XUP..TDOWN) as directions */
/* returning the coords of the 3rd nearest neighbor in that direction */

void third_neighbor(int x,int y,int z,int t,int *dirpt,int FB,int *xp,int *yp,int *zp,int *tp)
/* int x,y,z,t,*dirpt,FB;  coordinates of site, direction (eg XUP), and
				"forwards/backwards"  */
/* int *xp,*yp,*zp,*tp;    pointers to coordinates of neighbor */
{
  int dir;
  dir = ((FB == 1)? *dirpt : (7 -  *dirpt));
   *xp = x;
   *yp = y;
   *zp = z;
   *tp = t;
  switch(dir){
    case 0:
{
       *xp = ((x + 3) % nx);
      break; 
    }
    case 7:
{
       *xp = (((x + (4 * nx)) - 3) % nx);
      break; 
    }
    case 1:
{
       *yp = ((y + 3) % ny);
      break; 
    }
    case 6:
{
       *yp = (((y + (4 * ny)) - 3) % ny);
      break; 
    }
    case 2:
{
       *zp = ((z + 3) % nz);
      break; 
    }
    case 5:
{
       *zp = (((z + (4 * nz)) - 3) % nz);
      break; 
    }
    case 3:
{
       *tp = ((t + 3) % nt);
      break; 
    }
    case 4:
{
       *tp = (((t + (4 * nt)) - 3) % nt);
      break; 
    }
    default:
{
      printf("third_neighb: bad direction\n");
      exit(1);
    }
  }
}
