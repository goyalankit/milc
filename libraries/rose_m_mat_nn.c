/*******************  m_mat_nn.c  (in su3.a) ****************************
*									*
* void mult_su3_nn( su3_matrix *a,*b,*c )				*
* matrix multiply, no adjoints 						*
* C  <-  A*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_su3_nn(su3_matrix *a,su3_matrix *b,su3_matrix *c)
{
  register int i;
  register int j;
  register int k;
  register complex x;
  register complex y;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      x.real = (x.imag = 0.0);
//    __assume_aligned(&(a->e[0][0]), 64);
//    __assume_aligned(&(b->e[0][0]), 64);
//    __assume_aligned(&(c->e[0][0]), 64);
      
#pragma vector aligned always
      for (k = 0; k < 3; k++) {{
          y.real = (((a -> e)[i][k].real * (b -> e)[k][j].real) - ((a -> e)[i][k].imag * (b -> e)[k][j].imag));
          y.imag = (((a -> e)[i][k].real * (b -> e)[k][j].imag) + ((a -> e)[i][k].imag * (b -> e)[k][j].real));
        };
{
          x.real += y.real;
          x.imag += y.imag;
        };
      }
      (c -> e)[i][j] = x;
    }
  }
}
/* "Hand coded" routines, clearer coding is up above */
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */
#else
#endif	/* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#ifdef FAST" */
