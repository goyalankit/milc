/****************  m_mat_na.c  (in su3.a) *******************************
*									*
* void mult_su3_na( su3_matrix *a,*b,*c )				*
* matrix multiply, second matrix is adjoint 				*
* C  <-  A*B_adjoint							*
*/
#include <assert.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#ifndef FAST

void mult_su3_na(su3_matrix *a,su3_matrix *b,su3_matrix *c)
{
  int i;
  int j;
  int k;
  complex x;
  complex y;
  indigo__vector_stride_c(16,0,((void *)(&i)),sizeof(int ));
  indigo__vector_stride_c(16,0,((void *)(&i)),sizeof(int ));
  for (indigo__vector_stride_c(16,0,((void *)(&i)),sizeof(int )), i = 0; i < 3; i++) {
    indigo__vector_stride_c(16,1,((void *)(&j)),sizeof(int ));
    indigo__vector_stride_c(16,1,((void *)(&j)),sizeof(int ));
    for (indigo__vector_stride_c(16,1,((void *)(&j)),sizeof(int )), j = 0; j < 3; j++) {
      indigo__vector_stride_c(16,2,((void *)(&x.real)),sizeof(float ));
      indigo__vector_stride_c(16,3,((void *)(&x.imag)),sizeof(float ));
      x.real = (x.imag = 0.0);
//    __assume_aligned(&(a->e[0][0]), 64);
//    __assume_aligned(&(b->e[0][0]), 64);
//    __assume_aligned(&(c->e[0][0]), 64);
      
#pragma vector aligned always
      indigo__vector_stride_c(16,4,((void *)(&k)),sizeof(int ));
      indigo__vector_stride_c(16,4,((void *)(&k)),sizeof(int ));
      for (indigo__vector_stride_c(16,4,((void *)(&k)),sizeof(int )), k = 0; k < 3; k++) {{
          indigo__vector_stride_c(16,5,((void *)(&y.real)),sizeof(float ));
          indigo__vector_stride_c(16,6,((void *)(&(a -> e)[i][k].real)),sizeof(float ));
          indigo__vector_stride_c(16,7,((void *)(&(b -> e)[j][k].real)),sizeof(float ));
          indigo__vector_stride_c(16,8,((void *)(&(a -> e)[i][k].imag)),sizeof(float ));
          indigo__vector_stride_c(16,9,((void *)(&(b -> e)[j][k].imag)),sizeof(float ));
          y.real = (((a -> e)[i][k].real * (b -> e)[j][k].real) + ((a -> e)[i][k].imag * (b -> e)[j][k].imag));
          indigo__vector_stride_c(16,10,((void *)(&y.imag)),sizeof(float ));
          indigo__vector_stride_c(16,8,((void *)(&(a -> e)[i][k].imag)),sizeof(float ));
          indigo__vector_stride_c(16,7,((void *)(&(b -> e)[j][k].real)),sizeof(float ));
          indigo__vector_stride_c(16,6,((void *)(&(a -> e)[i][k].real)),sizeof(float ));
          indigo__vector_stride_c(16,9,((void *)(&(b -> e)[j][k].imag)),sizeof(float ));
          y.imag = (((a -> e)[i][k].imag * (b -> e)[j][k].real) - ((a -> e)[i][k].real * (b -> e)[j][k].imag));
        };
{
          indigo__vector_stride_c(16,2,((void *)(&x.real)),sizeof(float ));
          indigo__vector_stride_c(16,5,((void *)(&y.real)),sizeof(float ));
          x.real += y.real;
          indigo__vector_stride_c(16,3,((void *)(&x.imag)),sizeof(float ));
          indigo__vector_stride_c(16,10,((void *)(&y.imag)),sizeof(float ));
          x.imag += y.imag;
        };
      }
      indigo__vector_stride_c(16,11,((void *)(&(c -> e)[i][j])),sizeof(complex ));
      indigo__vector_stride_c(16,12,((void *)(&x)),sizeof(complex ));
      (c -> e)[i][j] = x;
    }
  }
}
/* "Hand coded" routines, clearer coding is up above */
#else
#endif	/* End of "#ifdef FAST" */

void indigo__create_map()
{
  indigo__write_idx_c("i",1);
  indigo__write_idx_c("j",1);
  indigo__write_idx_c("x.real",6);
  indigo__write_idx_c("x.imag",6);
  indigo__write_idx_c("k",1);
  indigo__write_idx_c("y.real",6);
  indigo__write_idx_c("(a -> e)[i][k].real",19);
  indigo__write_idx_c("(b -> e)[j][k].real",19);
  indigo__write_idx_c("(a -> e)[i][k].imag",19);
  indigo__write_idx_c("(b -> e)[j][k].imag",19);
  indigo__write_idx_c("y.imag",6);
  indigo__write_idx_c("(c -> e)[i][j]",14);
  indigo__write_idx_c("x",1);
}
