/*****************  sub4vecs.c  (in su3.a) ******************************
*									*
*  Subtract four su3_vectors from an su3_vector				*
* void sub_four_su3_vecs( su3_vector *a,*b1,*b2,*b3,*b4) 		*
* A  <-  A - B1 - B2 - B3 - B4						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* subtract four su3 vectors */
#ifndef FAST

void sub_four_su3_vecs(su3_vector *a,su3_vector *b1,su3_vector *b2,su3_vector *b3,su3_vector *b4)
{
  register int i;
  for (i = 0; i < 3; i++) {{
      (a -> c)[i].real = ((a -> c)[i].real - (b1 -> c)[i].real);
      (a -> c)[i].imag = ((a -> c)[i].imag - (b1 -> c)[i].imag);
    };
{
      (a -> c)[i].real = ((a -> c)[i].real - (b2 -> c)[i].real);
      (a -> c)[i].imag = ((a -> c)[i].imag - (b2 -> c)[i].imag);
    };
{
      (a -> c)[i].real = ((a -> c)[i].real - (b3 -> c)[i].real);
      (a -> c)[i].imag = ((a -> c)[i].imag - (b3 -> c)[i].imag);
    };
{
      (a -> c)[i].real = ((a -> c)[i].real - (b4 -> c)[i].real);
      (a -> c)[i].imag = ((a -> c)[i].imag - (b4 -> c)[i].imag);
    };
  }
}
#else
#endif
