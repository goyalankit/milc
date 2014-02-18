/*******************  trace_su3.c  (in su3.a) ***************************
*									*
* complex trace_su3(a) su3_matrix *a;					*
* return complex trace of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* Complex trace of an SU3 matrix */

complex trace_su3(su3_matrix *a)
{
  register complex t1;
  register complex t2;
{
    t1.real = ((a -> e)[0][0].real + (a -> e)[1][1].real);
    t1.imag = ((a -> e)[0][0].imag + (a -> e)[1][1].imag);
  };
{
    t2.real = (t1.real + (a -> e)[2][2].real);
    t2.imag = (t1.imag + (a -> e)[2][2].imag);
  };
  return t2;
}
