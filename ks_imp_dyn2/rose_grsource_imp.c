/************************ grsource_imp.c *****************************/
/* MIMD version 6 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif
/* construct a gaussian random vector, g_rand, and phi=M(dagger)*g_rand  */
/* "parity" is EVEN, ODD, or EVENANDODD.  The parity is the parity at
    which phi is computed.  g_rand must always be computed at all sites. */

void grsource_imp(field_offset dest,float mass,int parity)
{
  register int i;
  register int j;
  register site *s;
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    for (j = 0; j < 3; j++) {
#ifdef SITERAND
      s -> g_rand.c[j].real = gaussian_rand_no(&s -> site_prn);
      s -> g_rand.c[j].imag = gaussian_rand_no(&s -> site_prn);
#else
#endif
    }
  }
#ifdef FN
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
#endif
  dslash_fn(((field_offset )(((char *)(&lattice[0].g_rand)) - ((char *)(lattice + 0)))),dest,parity);
  scalar_mult_latvec(dest,(-1.0),dest,parity);
  scalar_mult_add_latvec(dest,((field_offset )(((char *)(&lattice[0].g_rand)) - ((char *)(lattice + 0)))),(2.0 * mass),dest,parity);
/* grsource */
}
/* Check congrad by multiplying src by M, compare result to g_rand */
/* Before calling checkmul() you should call grsource(EVENANDODD) and
   congrad(...,EVENANDODD) */

void checkmul_imp(field_offset src,float mass)
{
  register int i;
  register int j;
  register site *s;
  dslash_fn(src,((field_offset )(((char *)(&lattice[0].ttt)) - ((char *)(lattice + 0)))),3);
  scalar_mult_add_latvec(((field_offset )(((char *)(&lattice[0].ttt)) - ((char *)(lattice + 0)))),src,(2.0 * mass),((field_offset )(((char *)(&lattice[0].ttt)) - ((char *)(lattice + 0)))),3);
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    printf("Site %d %d %d %d\n",(s -> x),(s -> y),(s -> z),(s -> t));
    for (j = 0; j < 3; j++) {
      printf("%d %d\t%e\t%e\t%e\n",i,j,((double )s -> g_rand.c[j].real),((double )s -> ttt.c[j].real),(((double )s -> g_rand.c[j].real) - ((double )s -> ttt.c[j].real)));
      printf("%d %d\t%e\t%e\t%e\n",i,j,((double )s -> g_rand.c[j].imag),((double )s -> ttt.c[j].imag),(((double )s -> g_rand.c[j].imag) - ((double )s -> ttt.c[j].imag)));
    }
    printf("\n");
/**sleep(2);**/
  }
/* checkmul */
}
