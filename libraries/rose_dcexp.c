/********************** dcexp.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* double complex exponential */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dcexp(double_complex *a)
{
  double_complex c;
  double mag;
  mag = exp(((double )(a -> real)));
  c.real = (mag * cos(((double )(a -> imag))));
  c.imag = (mag * sin(((double )(a -> imag))));
  return c;
}
