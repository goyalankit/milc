/******************  s_m_mat.c  (in su3.a) ******************************
*									*
* void scalar_mult_su3_matrix( su3_matrix *a, float s, su3_matrix *b)	*
* B <- s*A								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* b <- s*a, matrices */

void scalar_mult_su3_matrix(su3_matrix *a,float s,su3_matrix *b)
{
#ifndef FAST
  register int i;
  register int j;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {
      (b -> e)[i][j].real = (s * (a -> e)[i][j].real);
      (b -> e)[i][j].imag = (s * (a -> e)[i][j].imag);
    }
#else
#ifdef NATIVEDOUBLE
#else
#endif
#endif
}
