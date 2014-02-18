/***************** m_amv_4vec.c in su3.a    *****************************
*									*
*  void mult_adj_su3_mat_4vec( su3_matrix *mat,			        *
*  su3_vector *src, su3_vector *dest0, *dest1, *dest2, *dest3 )		*
*  Multiply an su3_vector by an array of four adjoint su3_matrices,	*
*  result in four SEPARATE su3_vectors.		           		*
*  desti  <-  A_adjoint[i] * src					*
*  See also m_amv_4dir.c for the case desti = dest[i]                   *
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_adj_su3_mat_4vec(su3_matrix *mat,su3_vector *src,su3_vector *dest0,su3_vector *dest1,su3_vector *dest2,su3_vector *dest3)
{
  mult_adj_su3_mat_vec((mat + 0),src,dest0);
  mult_adj_su3_mat_vec((mat + 1),src,dest1);
  mult_adj_su3_mat_vec((mat + 2),src,dest2);
  mult_adj_su3_mat_vec((mat + 3),src,dest3);
}
#else
/* Fast code, with subroutines inlined */
#ifdef NATIVEDOUBLE
#else
#endif
#endif	/* End of "#ifndef FAST" */
