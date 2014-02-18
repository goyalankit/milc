/****************  s_m_sum_vec.c  (in su3.a) ****************************
*									*
* void scalar_mult_sum_su3_vector( su3_vector *a, su3_vector *b, float s )*
* A <- A + s*B,   A and B vectors 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* a <- a + s*b, vectors */

void scalar_mult_sum_su3_vector(su3_vector *a,su3_vector *b,float s)
{
#ifndef NATIVEDOUBLE
  register int i;
  for (i = 0; i < 3; i++) {
    (a -> c)[i].real += (s * (b -> c)[i].real);
    (a -> c)[i].imag += (s * (b -> c)[i].imag);
  }
#else /* RS6000 version */
#endif
}
