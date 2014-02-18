/********************  s_m_hwvec.c  (in su3.a) ********************
*
*void scalar_mult_hwvec(half_wilson_vector *src, float s,
	half_wilson_vector *dest)
*  Multiply a half Wilson vector by a scalar
* dest  <-  s*src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_hwvec(half_wilson_vector *src,float s,half_wilson_vector *dest)
{
#ifndef NATIVEDOUBLE
  register int i;
  for (i = 0; i < 2; i++) 
    scalar_mult_su3_vector(((src -> h) + i),s,((dest -> h) + i));
#else /* RS6000 version */
#endif
}
