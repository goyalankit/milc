/***************** su3_proj_w.c  (in su3.a) ****************************/
/* MIMD version 6 */
/*									*
* void su3_projector_w( wilson_vector *a, wilson_vector *b, su3_matrix *c )
* C  <- sum over spins of outer product of A.d[i] and B.d[i]		*
*  C_ij = sum( A_i * B_adjoint_j )					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST
//dune999
complex cmplx(float x,float y);

void su3_projector_w(wilson_vector *a,wilson_vector *b,su3_matrix *c)
{
  register int i;
  register int j;
  register int k;
  register complex cc;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {
      (c -> e)[i][j] = cmplx(0.0,0.0);
      for (k = 0; k < 4; k++) {{
          cc.real = (((a -> d)[k].c[i].real * (b -> d)[k].c[j].real) + ((a -> d)[k].c[i].imag * (b -> d)[k].c[j].imag));
          cc.imag = (((a -> d)[k].c[i].imag * (b -> d)[k].c[j].real) - ((a -> d)[k].c[i].real * (b -> d)[k].c[j].imag));
        };
{
          (c -> e)[i][j].real += cc.real;
          (c -> e)[i][j].imag += cc.imag;
        };
      }
    }
}
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */
#else
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
