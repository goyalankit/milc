/****************** ploop3.c ************************************/
/* MIMD version 6 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* DT 12/97 use local matrix "tmat" instead of "tempmat2" */
/* Uses only site structure member tempmat1 */
/* Macros ...
   BPCORR
   saves the Polyakov loop value for each site
   in the site member "ploop" for later use 
 */
#include "generic_includes.h"

complex ploop()
{
  register int i;
  register int t;
  register site *st;
  msg_tag *tag;
  complex sum;
  complex plp;
  su3_matrix tmat;
  int d[4UL];
  sum = cmplx(0.0,0.0);
  d[0] = (d[1] = (d[2] = 0));
/* First multiply the link on every even site by the link above it */
/* We will compute the Polyakov loop "at" the even sites in the 
	first two time slices. */
  tag = start_gather(((field_offset )(((char *)(lattice[0].link + 3)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),3,2,gen_pt[0]);
  wait_gather(tag);
  for (((i = 0) , (st = lattice)); i < even_sites_on_node; (i++ , st++)) {
    mult_su3_nn(((st -> link) + 3),((su3_matrix *)gen_pt[0][i]),&st -> tempmat1);
  }
  cleanup_gather(tag);
  for (t = 2; t < nt; t += 2) {
/* distance from which to gather */
    d[3] = t;
    tag = start_general_gather(((field_offset )(((char *)(&lattice[0].tempmat1)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),d,2,gen_pt[0]);
    wait_general_gather(tag);
    for (((i = 0) , (st = lattice)); i < even_sites_on_node; (i++ , st++)) {
/* only compute on first two slices */
      if ((st -> t) > 1) 
        continue; 
      mult_su3_nn(&st -> tempmat1,((su3_matrix *)gen_pt[0][i]),&tmat);
      lattice[i].tempmat1 = tmat;
/* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
    }
    cleanup_general_gather(tag);
  }
  for (((i = 0) , (st = lattice)); i < even_sites_on_node; (i++ , st++)) {
    if ((st -> t) > 1) 
      continue; 
    plp = trace_su3(&st -> tempmat1);
/* Save for later correlation measurements */
{
      sum.real += plp.real;
      sum.imag += plp.imag;
    };
#ifdef BPCORR
/* Save for subsequent correlation measurements */
/* Note the results are saved on even sites in
	   slices 0 and 1 */
#endif
  }
  g_complexsum(&sum);
  plp.real = (sum.real / ((float )((nx * ny) * nz)));
  plp.imag = (sum.imag / ((float )((nx * ny) * nz)));
  return plp;
}
