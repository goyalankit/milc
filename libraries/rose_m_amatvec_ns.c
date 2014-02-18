/******************  m_amatvec_ns.c  (in su3.a) *************************
*									*
* void mult_adj_su3_mat_vec_nsum( su3_matrix *a, su3_vector *b,*c )	*
* adjoint matrix times vector multiply and subtract from another vector *
* C  <-  C - A_adjoint*B						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_adj_su3_mat_vec_nsum(su3_matrix *a,su3_vector *b,su3_vector *c)
{
  register int i;
  register int j;
  register complex x;
  register complex y;
  register complex z;
  for (i = 0; i < 3; i++) {
    x.real = (x.imag = 0.0);
    for (j = 0; j < 3; j++) {{
        z.real = (a -> e)[j][i].real;
        z.imag = -(a -> e)[j][i].imag;
      };
{
        y.real = ((z.real * (b -> c)[j].real) - (z.imag * (b -> c)[j].imag));
        y.imag = ((z.real * (b -> c)[j].imag) + (z.imag * (b -> c)[j].real));
      }
{
        x.real += y.real;
        x.imag += y.imag;
      };
    }
    (c -> c)[i].real -= x.real;
    (c -> c)[i].imag -= x.imag;
  }
}
#else
#ifdef NATIVEDOUBLE
#else
#endif
#endif	/* End of "#ifdef FAST" */
