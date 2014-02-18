/**************  m_su2_mat_vec_a.c (in su3.a) **********************
*									*
*  adjoint su2 matrix times vector                             		*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1)
{
/* Multiplies the complex row spinor (x0, x1) by the adjoint of the */
/* SU(2) matrix u and puts the result in (x0,x1).  */
/* Thus x <-  x * u-adj       */
/* C. DeTar 3 Oct 1990 */
  complex z0;
  complex z1;
  complex t0;
  complex t1;
  t0 =  *x0;
  t1 =  *x1;
{
    z0.real = ((t0.real * (u -> e)[0][0].real) + (t0.imag * (u -> e)[0][0].imag));
    z0.imag = ((t0.imag * (u -> e)[0][0].real) - (t0.real * (u -> e)[0][0].imag));
  };
{
    z1.real = ((t1.real * (u -> e)[0][1].real) + (t1.imag * (u -> e)[0][1].imag));
    z1.imag = ((t1.imag * (u -> e)[0][1].real) - (t1.real * (u -> e)[0][1].imag));
  };
{
    x0 -> real = (z0.real + z1.real);
    x0 -> imag = (z0.imag + z1.imag);
  };
{
    z0.real = ((t0.real * (u -> e)[1][0].real) + (t0.imag * (u -> e)[1][0].imag));
    z0.imag = ((t0.imag * (u -> e)[1][0].real) - (t0.real * (u -> e)[1][0].imag));
  };
{
    z1.real = ((t1.real * (u -> e)[1][1].real) + (t1.imag * (u -> e)[1][1].imag));
    z1.imag = ((t1.imag * (u -> e)[1][1].real) - (t1.real * (u -> e)[1][1].imag));
  };
{
    x1 -> real = (z0.real + z1.real);
    x1 -> imag = (z0.imag + z1.imag);
  };
/* m_su2_mat_vec_a.c */
}
