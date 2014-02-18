/*************************** ranmom.c *******************************/
/* MIMD version 6 */
/* Produce Gaussian random momenta for the gauge fields. */
#include "generic_includes.h"
#include <defines.h>                 /* For SITERAND */

void ranmom()
{
  register int i;
  register int dir;
  register site *s;
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    for (dir = 0; dir <= 3; dir++) {
#ifdef SCHROED_FUN
#endif
#ifdef SITERAND
      random_anti_hermitian(((s -> mom) + dir),&s -> site_prn);
#else
#endif
#ifdef SCHROED_FUN
#endif
    }
  }
}
