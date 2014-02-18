/*****************  m_amv_4dir.c  (in su3.a) *****************************
*									*
*  void mult_adj_su3_mat_vec_4dir( su3_matrix *mat,			*
*  su3_vector *src, su3_vector *dest )					*
*  Multiply an su3_vector by an array of four adjoint su3_matrices,	*
*  result in an array of four su3_vectors.				*
*  dest[i]  <-  A_adjoint[i] * src					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_adj_su3_mat_vec_4dir(su3_matrix *mat,su3_vector *src,su3_vector *dest)
{
  mult_adj_su3_mat_vec((mat + 0),src,(dest + 0));
  mult_adj_su3_mat_vec((mat + 1),src,(dest + 1));
  mult_adj_su3_mat_vec((mat + 2),src,(dest + 2));
  mult_adj_su3_mat_vec((mat + 3),src,(dest + 3));
}
#else
/* Fast code, with subroutines inlined */
#ifdef NATIVEDOUBLE
#else
#endif
#endif	/* End of "#ifndef FAST" */
