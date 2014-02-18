/****************  cs_m_s_mat.c  (in su3.a) *****************************
*									*
* void c_scalar_mult_sub_su3mat( su3_matrix *a, su3_matrix *b,		*
*	complex *s, su3_matrix *c)					*
* C <- A - s*B,   A,B and C matrices 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* c <- a - s*b, matrices */

void c_scalar_mult_sub_su3mat(su3_matrix *a,su3_matrix *b,complex *s,su3_matrix *c)
{
#ifndef NATIVEDOUBLE
  register int i;
  register int j;
  complex t;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {
      t = cmul(((b -> e)[i] + j),s);
      (c -> e)[i][j] = csub(((a -> e)[i] + j),&t);
    }
#else
#endif
}
