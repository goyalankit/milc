/*****************  wvec_rdot.c  (in su3.a) ******************************
*									*
* float wvec_rdot( wilson_vector *a, wilson_vector *b )			*
* return real part of dot product of two wilson_vectors			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

float wvec_rdot(wilson_vector *a,wilson_vector *b)
{
#ifndef FAST
  register float temp1;
  register float temp2;
  register int i;
  temp2 = 0.0;
  for (i = 0; i < 4; i++) {
    temp1 = ((a -> d)[i].c[0].real * (b -> d)[i].c[0].real);
    temp2 += temp1;
    temp1 = ((a -> d)[i].c[0].imag * (b -> d)[i].c[0].imag);
    temp2 += temp1;
    temp1 = ((a -> d)[i].c[1].real * (b -> d)[i].c[1].real);
    temp2 += temp1;
    temp1 = ((a -> d)[i].c[1].imag * (b -> d)[i].c[1].imag);
    temp2 += temp1;
    temp1 = ((a -> d)[i].c[2].real * (b -> d)[i].c[2].real);
    temp2 += temp1;
    temp1 = ((a -> d)[i].c[2].imag * (b -> d)[i].c[2].imag);
    temp2 += temp1;
  }
  return temp2;
#else
#ifndef NATIVEDOUBLE
#else
#endif
#endif
}
