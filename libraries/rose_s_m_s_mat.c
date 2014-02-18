/****************  s_m_s_mat.c  (in su3.a) ******************************
*									*
* void scalar_mult_sub_su3_matrix( su3_matrix *a, su3_matrix *b,	*
*	float s, su3_matrix *c)						*
* C <- A - s*B,   A,B and C matrices 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* c <- a - s*b, matrices */

void scalar_mult_sub_su3_matrix(su3_matrix *a,su3_matrix *b,float s,su3_matrix *c)
{
#ifndef NATIVEDOUBLE
  register int i;
  register int j;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {
      (c -> e)[i][j].real = ((a -> e)[i][j].real - (s * (b -> e)[i][j].real));
      (c -> e)[i][j].imag = ((a -> e)[i][j].imag - (s * (b -> e)[i][j].imag));
    }
#else /* RS6000 version */
#endif
}
