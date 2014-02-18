/********************  msq_wvec.c  (in su3.a) ********************
*
*float msq_wvec(wilson_vector *vec)
*  squared magnitude of a Wilson vector
* 
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

float magsq_wvec(wilson_vector *vec)
{
  register int i;
  register float sum;
  sum = 0.0;
  for (i = 0; i < 4; i++) 
    sum += magsq_su3vec(((vec -> d) + i));
  return sum;
#else /* Fast version */
#ifdef NATIVEDOUBLE
#else
#endif
#endif
}
