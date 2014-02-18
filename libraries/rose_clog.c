/********************** clog.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* complex logarithm */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex clog(complex *a)
{
  complex c;
  c.real = (0.5 * ((float )(log(((double )(((a -> real) * (a -> real)) + ((a -> imag) * (a -> imag))))))));
  c.imag = ((float )(atan2(((double )(a -> imag)),((double )(a -> real)))));
  return c;
}
