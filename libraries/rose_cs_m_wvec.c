/********************  	cs_m_wvec.c  (in su3.a) ********************
*
*void c_scalar_mult_wvec(wilson_vector *src, complex *s, wilson_vector *dest)
*  Multiply a Wilson vector by a complex scalar and add to another vector
* dest  <-  s * src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_wvec(wilson_vector *src,complex *phase,wilson_vector *dest)
{
#ifndef FAST
  register int i;
  register int j;
  complex t;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {{
        (dest -> d)[i].c[j].real = (((src -> d)[i].c[j].real * (phase -> real)) - ((src -> d)[i].c[j].imag * (phase -> imag)));
        (dest -> d)[i].c[j].imag = (((src -> d)[i].c[j].real * (phase -> imag)) + ((src -> d)[i].c[j].imag * (phase -> real)));
      };
    }
  }
#else
#ifdef NATIVEDOUBLE
#else
#endif
#endif
}
