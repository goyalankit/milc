/****************  m_matvec.c  (in su3.a) *******************************
*									*
* void mult_su3_mat_vec( su3_matrix *a, su3_vector *b,*c )		*
* matrix times vector multiply, no adjoints 				*
*  C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_su3_mat_vec(su3_matrix *a,su3_vector *b,su3_vector *c)
{
  register int i;
  register int j;
  register complex x;
  register complex y;
  for (i = 0; i < 3; i++) {
    x.real = (x.imag = 0.0);
//    __assume_aligned(&(a->e[0][0]), 64);
//    __assume_aligned(&(b->c[0]), 64);
//    __assume_aligned(&(c->c[0]), 64);
    
#pragma vector always aligned
    for (j = 0; j < 3; j++) {{
        y.real = (((a -> e)[i][j].real * (b -> c)[j].real) - ((a -> e)[i][j].imag * (b -> c)[j].imag));
        y.imag = (((a -> e)[i][j].real * (b -> c)[j].imag) + ((a -> e)[i][j].imag * (b -> c)[j].real));
      }
{
        x.real += y.real;
        x.imag += y.imag;
      };
    }
    (c -> c)[i] = x;
  }
}
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */
#else
#endif	/* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#infdef FAST" */
