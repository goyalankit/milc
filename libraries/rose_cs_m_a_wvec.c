/********************  	cs_m_a_wvec.c  (in su3.a) ********************
*
*void c_scalar_mult_add_wvec(wilson_vector *src1, wilson_vector *src2,
	complex *s, wilson_vector *dest)
*  Multiply a Wilson vector by a complex scalar and add to another vector
* dest  <-  src1 + s*src2
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_wvec(wilson_vector *src1,wilson_vector *src2,complex *phase,wilson_vector *dest)
{
  register int i;
  register int j;
  register float sr;
  register float si;
  register float br;
  register float bi;
  register float cr;
  register float ci;
  sr = (phase -> real);
  si = (phase -> imag);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 3; j++) {
      br = (src2 -> d)[i].c[j].real;
      bi = (src2 -> d)[i].c[j].imag;
      cr = ((sr * br) - (si * bi));
      ci = ((sr * bi) + (si * br));
      (dest -> d)[i].c[j].real = ((src1 -> d)[i].c[j].real + cr);
      (dest -> d)[i].c[j].imag = ((src1 -> d)[i].c[j].imag + ci);
    }
  }
}
