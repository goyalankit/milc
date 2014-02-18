/*****************  m_amatvec.c  (in su3.a) *****************************
*									*
*  void mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b,*c )		*
*  C  <-  A_adjoint * B							*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST
/* adjoint matrix times vector multiply */

void mult_adj_su3_mat_vec(su3_matrix *a,su3_vector *b,su3_vector *c)
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
    (c -> c)[i] = x;
  }
}
#else
#ifdef NATIVEDOUBLE /* IBM RS6000 version */
#else
#endif	/* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#ifndef FAST" */
