/*****************  su3_rdot.c  (in su3.a) ******************************
*									*
* float su3_rdot( su3_vector *a, su3_vector *b )			*
* return real part of dot product of two su3_vectors			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

float su3_rdot(su3_vector *a,su3_vector *b)
{
#ifndef NATIVEDOUBLE
  register float temp1;
  register float temp2;
  temp2 = ((a -> c)[0].real * (b -> c)[0].real);
  temp1 = ((a -> c)[0].imag * (b -> c)[0].imag);
  temp2 += temp1;
  temp1 = ((a -> c)[1].real * (b -> c)[1].real);
  temp2 += temp1;
  temp1 = ((a -> c)[1].imag * (b -> c)[1].imag);
  temp2 += temp1;
  temp1 = ((a -> c)[2].real * (b -> c)[2].real);
  temp2 += temp1;
  temp1 = ((a -> c)[2].imag * (b -> c)[2].imag);
  temp2 += temp1;
  return temp2;
#else /* RS6000 version */
#endif
}
