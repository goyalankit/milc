/*****************  make_ahmat.c  (in su3.a) ****************************
*									*
* void make_anti_hermitian( su3_matrix *m3, anti_hermitmat *ah3)	*
* take the traceless and anti_hermitian part of an su3 matrix 		*
* and compress it 							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void make_anti_hermitian(su3_matrix *m3,anti_hermitmat *ah3)
{
  float temp;
  temp = ((((m3 -> e)[0][0].imag + (m3 -> e)[1][1].imag) + (m3 -> e)[2][2].imag) * 0.33333333);
  ah3 -> m00im = ((m3 -> e)[0][0].imag - temp);
  ah3 -> m11im = ((m3 -> e)[1][1].imag - temp);
  ah3 -> m22im = ((m3 -> e)[2][2].imag - temp);
  ah3 -> m01.real = (((m3 -> e)[0][1].real - (m3 -> e)[1][0].real) * 0.5);
  ah3 -> m02.real = (((m3 -> e)[0][2].real - (m3 -> e)[2][0].real) * 0.5);
  ah3 -> m12.real = (((m3 -> e)[1][2].real - (m3 -> e)[2][1].real) * 0.5);
  ah3 -> m01.imag = (((m3 -> e)[0][1].imag + (m3 -> e)[1][0].imag) * 0.5);
  ah3 -> m02.imag = (((m3 -> e)[0][2].imag + (m3 -> e)[2][0].imag) * 0.5);
  ah3 -> m12.imag = (((m3 -> e)[1][2].imag + (m3 -> e)[2][1].imag) * 0.5);
/* make_anti_hermitian */
}
#else
/* make_anti_hermitian */
#endif /*end ifdef FAST */