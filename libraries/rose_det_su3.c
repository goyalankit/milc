/******************  det_su3.c  (in su3.a) ******************************
*									*
* complex det_su3( su3_matrix *a )					*
* Complex determinant of an SU3 matrix 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
/* FIX THIS - more efficient to take cross product of first two
   rows, dot with third. */

complex det_su3(su3_matrix *a)
{
  register complex cc;
  register complex dd;
  register complex sum;
{
    cc.real = (((a -> e)[0][0].real * (a -> e)[1][1].real) - ((a -> e)[0][0].imag * (a -> e)[1][1].imag));
    cc.imag = (((a -> e)[0][0].real * (a -> e)[1][1].imag) + ((a -> e)[0][0].imag * (a -> e)[1][1].real));
  };
{
    sum.real = ((cc.real * (a -> e)[2][2].real) - (cc.imag * (a -> e)[2][2].imag));
    sum.imag = ((cc.real * (a -> e)[2][2].imag) + (cc.imag * (a -> e)[2][2].real));
  };
{
    cc.real = (((a -> e)[0][0].real * (a -> e)[1][2].real) - ((a -> e)[0][0].imag * (a -> e)[1][2].imag));
    cc.imag = (((a -> e)[0][0].real * (a -> e)[1][2].imag) + ((a -> e)[0][0].imag * (a -> e)[1][2].real));
  };
{
    dd.real = ((cc.real * (a -> e)[2][1].real) - (cc.imag * (a -> e)[2][1].imag));
    dd.imag = ((cc.real * (a -> e)[2][1].imag) + (cc.imag * (a -> e)[2][1].real));
  };
{
    sum.real = (sum.real - dd.real);
    sum.imag = (sum.imag - dd.imag);
  };
{
    cc.real = (((a -> e)[0][1].real * (a -> e)[1][2].real) - ((a -> e)[0][1].imag * (a -> e)[1][2].imag));
    cc.imag = (((a -> e)[0][1].real * (a -> e)[1][2].imag) + ((a -> e)[0][1].imag * (a -> e)[1][2].real));
  };
{
    dd.real = ((cc.real * (a -> e)[2][0].real) - (cc.imag * (a -> e)[2][0].imag));
    dd.imag = ((cc.real * (a -> e)[2][0].imag) + (cc.imag * (a -> e)[2][0].real));
  };
{
    sum.real = (sum.real + dd.real);
    sum.imag = (sum.imag + dd.imag);
  };
{
    cc.real = (((a -> e)[0][1].real * (a -> e)[1][0].real) - ((a -> e)[0][1].imag * (a -> e)[1][0].imag));
    cc.imag = (((a -> e)[0][1].real * (a -> e)[1][0].imag) + ((a -> e)[0][1].imag * (a -> e)[1][0].real));
  };
{
    dd.real = ((cc.real * (a -> e)[2][2].real) - (cc.imag * (a -> e)[2][2].imag));
    dd.imag = ((cc.real * (a -> e)[2][2].imag) + (cc.imag * (a -> e)[2][2].real));
  };
{
    sum.real = (sum.real - dd.real);
    sum.imag = (sum.imag - dd.imag);
  };
{
    cc.real = (((a -> e)[0][2].real * (a -> e)[1][0].real) - ((a -> e)[0][2].imag * (a -> e)[1][0].imag));
    cc.imag = (((a -> e)[0][2].real * (a -> e)[1][0].imag) + ((a -> e)[0][2].imag * (a -> e)[1][0].real));
  };
{
    dd.real = ((cc.real * (a -> e)[2][1].real) - (cc.imag * (a -> e)[2][1].imag));
    dd.imag = ((cc.real * (a -> e)[2][1].imag) + (cc.imag * (a -> e)[2][1].real));
  };
{
    sum.real = (sum.real + dd.real);
    sum.imag = (sum.imag + dd.imag);
  };
{
    cc.real = (((a -> e)[0][2].real * (a -> e)[1][1].real) - ((a -> e)[0][2].imag * (a -> e)[1][1].imag));
    cc.imag = (((a -> e)[0][2].real * (a -> e)[1][1].imag) + ((a -> e)[0][2].imag * (a -> e)[1][1].real));
  };
{
    dd.real = ((cc.real * (a -> e)[2][0].real) - (cc.imag * (a -> e)[2][0].imag));
    dd.imag = ((cc.real * (a -> e)[2][0].imag) + (cc.imag * (a -> e)[2][0].real));
  };
{
    sum.real = (sum.real - dd.real);
    sum.imag = (sum.imag - dd.imag);
  };
  return sum;
}
