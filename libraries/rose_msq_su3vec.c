/****************** msq_su3vec.c  (in su3.a) ****************************/
/* MIMD version 6 */
/*									*
* float magsq_su3vec( su3_vector *a )					*
* return squared magnitude of an SU3 vector
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

float magsq_su3vec(su3_vector *a)
{
  register float sum;
  register int i;
  for (((i = 0) , (sum = 0.0)); i < 3; i++) 
    sum += (((a -> c)[i].real * (a -> c)[i].real) + ((a -> c)[i].imag * (a -> c)[i].imag));
  return sum;
}
#else
#ifdef NATIVEDOUBLE /* IBM RS6000 version */
#else
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
