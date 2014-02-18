/**************** f_meas.c ***************************************/
/* MIMD version 6 */
/* UH 11/01 Include measurements for quark number susceptibilities */
/*          Note: this does not work for p4-action! */
/* UH 11/1/01 write complex stochastic estimators */
/* TB 10/01 Include measurements of dM/du0 for EOS */
/* CD 7/14/01 allow for multiple stochastic estimators NPBP_REPS */
/* DT 12/97 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
/* Measure fermionic observables:
    psi-bar-psi (separately on even and odd sites)
    fermion action
    This routine uses g_rand, phi, xxx, and other vectors used by
    the matrix inversion.
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */

void f_meas_imp(field_offset phi_off,field_offset xxx_off,float mass)
{
  float r_psi_bar_psi_even;
  float i_psi_bar_psi_even;
  float r_psi_bar_psi_odd;
  float i_psi_bar_psi_odd;
  float r_ferm_action;
/* local variables for accumulators */
  register int i;
  register site *st;
  double rfaction;
  double_complex pbp_e;
  double_complex pbp_o;
  double pbp_pbp;
  complex cc;
#ifdef DM_DU0
#endif
#ifdef CHEM_POT
#ifndef FN		/* FN is assumed for quark number susc. */
#endif
#ifndef NPBP_REPS	/* Need multiple repetitions for susceptibilities! */
#endif
#endif
/* If this feature is used more commonly, we should make npbp_reps
       a user-supplied parameter, instead of a macro and define it
       globally */
#ifdef NPBP_REPS
/* Number of repetitions of stochastic
                                   estimate */
#else
  int npbp_reps = 1;
#endif
  int jpbp_reps;
#ifdef FN
  if (!(valid_longlinks == 1)) 
    load_longlinks();
  if (!(valid_fatlinks == 1)) 
    load_fatlinks();
#endif
  for (jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++) {
    rfaction = ((double )0.0);
    pbp_e = (pbp_o = dcmplx(((double )0.0),((double )0.0)));
/* Make random source, and do inversion */
    grsource_imp(phi_off,mass,3);
    mat_invert_uml(((field_offset )(((char *)(&lattice[0].g_rand)) - ((char *)(lattice + 0)))),xxx_off,phi_off,mass);
#ifdef DM_DU0
#endif
#ifdef CHEM_POT
/* Start gathers from positive t-direction */
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* Start gathers from negative t-direction */
/* Wait gathers from positive t-direction and multiply by matrix */
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* Wait gathers from negative t-direction */
#endif
/* fermion action = phi.xxx */
/* psi-bar-psi on even sites = g_rand.xxx */
    for (((i = 0) , (st = lattice)); i < even_sites_on_node; (i++ , st++)) {
      cc = su3_dot(((su3_vector *)(((char *)st) + phi_off)),((su3_vector *)(((char *)st) + xxx_off)));
      rfaction += cc.real;
      cc = su3_dot(&st -> g_rand,((su3_vector *)(((char *)st) + xxx_off)));
{
        pbp_e.real += cc.real;
        pbp_e.imag += cc.imag;
      };
#ifdef DM_DU0
#endif
#ifdef CHEM_POT
#endif
    }
/* psi-bar-psi on odd sites */
    for (((i = even_sites_on_node) , (st = (lattice + i))); i < sites_on_node; (i++ , st++)) {
      cc = su3_dot(&st -> g_rand,((su3_vector *)(((char *)st) + xxx_off)));
{
        pbp_o.real += cc.real;
        pbp_o.imag += cc.imag;
      };
#ifdef DM_DU0
#endif
#ifdef CHEM_POT
#endif
    }
    g_dcomplexsum(&pbp_o);
    g_dcomplexsum(&pbp_e);
    g_doublesum(&rfaction);
#ifdef DM_DU0
#endif
    r_psi_bar_psi_odd = (pbp_o.real * (2.0 / ((double )volume)));
    i_psi_bar_psi_odd = (pbp_o.imag * (2.0 / ((double )volume)));
    r_psi_bar_psi_even = (pbp_e.real * (2.0 / ((double )volume)));
    i_psi_bar_psi_even = (pbp_e.imag * (2.0 / ((double )volume)));
    r_ferm_action = (rfaction * (1.0 / ((double )volume)));
    if (this_node == 0) 
      printf("PBP: mass %e     %e  %e  %e  %e ( %d of %d )\n",mass,r_psi_bar_psi_even,r_psi_bar_psi_odd,i_psi_bar_psi_even,i_psi_bar_psi_odd,(jpbp_reps + 1),npbp_reps);
    if (this_node == 0) 
      printf("FACTION: mass = %e,  %e ( %d of %d )\n",mass,r_ferm_action,(jpbp_reps + 1),npbp_reps);
#ifdef CHEM_POT
/* free up the buffers */
#endif
#ifdef NPBP_REPS
#endif
#ifdef CHEM_POT
/* Start gathers from positive t-direction */
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* Start gathers from negative t-direction */
/* Wait gathers from positive t-direction and multiply by matrix */
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* Wait gathers from negative t-direction */
/* free up the buffers */
#endif
  }
}
