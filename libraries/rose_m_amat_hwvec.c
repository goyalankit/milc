/**************  m_amat_hwvec.c  (in su3.a) **********************
*									*
*  void mult_adj_su3_mat_hwvec( su3_matrix *mat,			*
*	half_wilson_vector *src,*dest )					*
*  multiply a Wilson half-vector by the adjoint of a matrix		*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_adj_su3_mat_hwvec(su3_matrix *mat,half_wilson_vector *src,half_wilson_vector *dest)
{
  mult_adj_su3_mat_vec(mat,((src -> h) + 0),((dest -> h) + 0));
  mult_adj_su3_mat_vec(mat,((src -> h) + 1),((dest -> h) + 1));
}
#else  /* Fast version */
#ifdef NATIVEDOUBLE
#else
#endif
/*    mult_adj_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */
/*    mult_adj_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */
#endif /* "ifndef FAST" */
