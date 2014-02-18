/********************  addvec.c  (in su3.a) *****************************
*									*
*  Add two SU3 vectors							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_vector(su3_vector *a,su3_vector *b,su3_vector *c)
{
  register int i;
  for (i = 0; i < 3; i++) {{
      (c -> c)[i].real = ((a -> c)[i].real + (b -> c)[i].real);
      (c -> c)[i].imag = ((a -> c)[i].imag + (b -> c)[i].imag);
    };
  }
}
