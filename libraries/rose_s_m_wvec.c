/********************  s_m_wvec.c  (in su3.a) ********************
*
*void scalar_mult_wvec(wilson_vector *src, float s, wilson_vector *dest)
*  Multiply a Wilson vector by a scalar
* dest  <-  s*src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_wvec(wilson_vector *src,float s,wilson_vector *dest)
{
#ifndef NATIVEDOUBLE
  register int i;
  for (i = 0; i < 4; i++) 
    scalar_mult_su3_vector(((src -> d) + i),s,((dest -> d) + i));
#else /* RS6000 version */
#endif
}
