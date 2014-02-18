/********************* rephase.c *******************************/
#include "generic_ks_includes.h"
/********* phaseset() - set up KS phase vectors **********/
/* ANTIPERIODIC bc's in t and PERIODIC in x,y,z */

void phaseset()
{
/* pointer to current site */
  register site *sit;
  register int i;
/*	phase of link(i,mu) = sum(nu<mu) { -1^i[nu] }		*/
/*	all t phases for t=nt-1 time slice get extra minus sign	*/
/*	   to give antiperiodic boundary conditions		*/
  for (((i = 0) , (sit = lattice)); i < sites_on_node; (i++ , sit++)) {
    (sit -> phase)[3] = 1.0;
    if (((sit -> t) % 2) == 1) 
      (sit -> phase)[0] = (-1.0);
    else 
      (sit -> phase)[0] = 1.0;
    if (((sit -> x) % 2) == 1) 
      (sit -> phase)[1] = -(sit -> phase)[0];
    else 
      (sit -> phase)[1] = (sit -> phase)[0];
    if (((sit -> y) % 2) == 1) 
      (sit -> phase)[2] = -(sit -> phase)[1];
    else 
      (sit -> phase)[2] = (sit -> phase)[1];
    if ((sit -> t) == (nt - 1)) {
/* antiperiodic boundary conditions in Euclidean time */
      (sit -> phase)[3] = -(sit -> phase)[3];
    }
  }
}
/************************** rephase() ******************************/
/* put Kogut-Sussind and boundary condition phase factors into or
   out of lattice */

void rephase(int flag)
{
  register int i;
  register int j;
  register int k;
  register int dir;
  register site *s;
/* Check to make sure we are going in expected direction */
  if (!(((flag == 1) && (phases_in == 0)) || ((flag == 0) && (phases_in == 1)))) {
    if (this_node == 0) 
      printf("DUMMY: you fouled up the phases\n");
    terminate(1);
  }
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    for (dir = 0; dir <= 3; dir++) {
      for (j = 0; j < 3; j++) 
        for (k = 0; k < 3; k++) {
          (s -> link)[dir].e[j][k].real *= (s -> phase)[dir];
          (s -> link)[dir].e[j][k].imag *= (s -> phase)[dir];
        }
    }
  }
  phases_in = flag;
/* rephase.c */
}
