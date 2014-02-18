/*****************  s_m_atm_wvec.c  (in su3.a) ********************
*
*void scalar_mult_addtm_wvec(wilson_vector *src1, wilson_vector *src2,
	float s, wilson_vector *dest)
*  Multiply a Wilson vector by a scalar and add to minus one times
*   another vector
* dest  <-  (-1)*src1 + s*src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_addtm_wvec(wilson_vector *src1,wilson_vector *src2,float s,wilson_vector *dest)
{
#ifndef NATIVEDOUBLE
  register int i;
  register int j;
/*spins*/
  for (i = 0; i < 4; i++) {
/*colors*/
    for (j = 0; j < 3; j++) {
      (dest -> d)[i].c[j].real = (-(src1 -> d)[i].c[j].real + (s * (src2 -> d)[i].c[j].real));
      (dest -> d)[i].c[j].imag = (-(src1 -> d)[i].c[j].imag + (s * (src2 -> d)[i].c[j].imag));
    }
  }
#else /* RS6000 version */
#endif
}
