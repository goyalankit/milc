/******************  wvec_dot.c  (in su3.a) ****************************/
/* MIMD version 6 */
/*									*
* complex wvec_dot(a,b) wilson_vector *a,*b;				*
* return dot product of two wilson_vectors					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex wvec_dot(wilson_vector *a,wilson_vector *b)
{
#ifndef FAST
  complex temp1;
  complex temp2;
  register int i;
  temp1.real = (temp1.imag = 0.0);
  for (i = 0; i < 4; i++) {{
      temp2.real = (((a -> d)[i].c[0].real * (b -> d)[i].c[0].real) + ((a -> d)[i].c[0].imag * (b -> d)[i].c[0].imag));
      temp2.imag = (((a -> d)[i].c[0].real * (b -> d)[i].c[0].imag) - ((a -> d)[i].c[0].imag * (b -> d)[i].c[0].real));
    };
{
      temp1.real += temp2.real;
      temp1.imag += temp2.imag;
    };
{
      temp2.real = (((a -> d)[i].c[1].real * (b -> d)[i].c[1].real) + ((a -> d)[i].c[1].imag * (b -> d)[i].c[1].imag));
      temp2.imag = (((a -> d)[i].c[1].real * (b -> d)[i].c[1].imag) - ((a -> d)[i].c[1].imag * (b -> d)[i].c[1].real));
    };
{
      temp1.real += temp2.real;
      temp1.imag += temp2.imag;
    };
{
      temp2.real = (((a -> d)[i].c[2].real * (b -> d)[i].c[2].real) + ((a -> d)[i].c[2].imag * (b -> d)[i].c[2].imag));
      temp2.imag = (((a -> d)[i].c[2].real * (b -> d)[i].c[2].imag) - ((a -> d)[i].c[2].imag * (b -> d)[i].c[2].real));
    };
{
      temp1.real += temp2.real;
      temp1.imag += temp2.imag;
    };
  }
  return temp1;
#else
#ifdef NATIVEDOUBLE
#else
#endif
#endif
}
