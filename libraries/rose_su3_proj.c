/*****************  su3_proj.c  (in su3.a) ******************************
*									*
* void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c )	*
* C  <- outer product of A and B					*
*  C_ij = A_i * B_adjoint_j						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void su3_projector(su3_vector *a,su3_vector *b,su3_matrix *c)
{
  register int i;
  register int j;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {{
        (c -> e)[i][j].real = (((a -> c)[i].real * (b -> c)[j].real) + ((a -> c)[i].imag * (b -> c)[j].imag));
        (c -> e)[i][j].imag = (((a -> c)[i].imag * (b -> c)[j].real) - ((a -> c)[i].real * (b -> c)[j].imag));
      };
    }
}
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */
#else
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif /* end ifdef FAST */
