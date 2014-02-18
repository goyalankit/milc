/********************** ce_itheta.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* exp( i*theta ) */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex ce_itheta(float theta)
{
  complex c;
  c.real = ((float )(cos(((double )theta))));
  c.imag = ((float )(sin(((double )theta))));
/* there must be a more efficient way */
  return c;
}
