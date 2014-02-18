/******************  su3_dot.c  (in su3.a) ******************************
*									*
* complex su3_dot( su3_vector *a, su3_vector *b )			*
* return dot product of two su3_vectors: a^dagger b			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex su3_dot(su3_vector *a,su3_vector *b)
{
#ifndef FAST
  complex temp1;
  complex temp2;
{
    temp1.real = (((a -> c)[0].real * (b -> c)[0].real) + ((a -> c)[0].imag * (b -> c)[0].imag));
    temp1.imag = (((a -> c)[0].real * (b -> c)[0].imag) - ((a -> c)[0].imag * (b -> c)[0].real));
  }
{
    temp2.real = (((a -> c)[1].real * (b -> c)[1].real) + ((a -> c)[1].imag * (b -> c)[1].imag));
    temp2.imag = (((a -> c)[1].real * (b -> c)[1].imag) - ((a -> c)[1].imag * (b -> c)[1].real));
  }
{
    temp1.real += temp2.real;
    temp1.imag += temp2.imag;
  };
{
    temp2.real = (((a -> c)[2].real * (b -> c)[2].real) + ((a -> c)[2].imag * (b -> c)[2].imag));
    temp2.imag = (((a -> c)[2].real * (b -> c)[2].imag) - ((a -> c)[2].imag * (b -> c)[2].real));
  }
{
    temp1.real += temp2.real;
    temp1.imag += temp2.imag;
  };
  return temp1;
#else /* RS6000 version */
#ifdef NATIVEDOUBLE
#else
#endif
#endif
}
