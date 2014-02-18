/****************  m_mv_s_4dir.c  (in su3.a) *****************************
*									*
* void mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b[0123],*c )*
* Multiply the elements of an array of four su3_matrices by the		*
* four su3_vectors, and add the results to				*
* produce a single su3_vector.						*
* C  <-  A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3]			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_su3_mat_vec_sum_4dir(su3_matrix *a,su3_vector *b0,su3_vector *b1,su3_vector *b2,su3_vector *b3,su3_vector *c)
{
  mult_su3_mat_vec((a + 0),b0,c);
  mult_su3_mat_vec_sum((a + 1),b1,c);
  mult_su3_mat_vec_sum((a + 2),b2,c);
  mult_su3_mat_vec_sum((a + 3),b3,c);
}
#else
/* Fast code, with subroutines inlined */
#ifdef NATIVEDOUBLE /* IBM RS6000 version */
#else
#endif  /* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#ifdef FAST" */
