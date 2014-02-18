/********************** dce_itheta.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* double complex exp( i*theta ) */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

double_complex dce_itheta(double theta)
{
  double_complex c;
  c.real = cos(((double )theta));
  c.imag = sin(((double )theta));
/* there must be a more efficient way */
  return c;
}
