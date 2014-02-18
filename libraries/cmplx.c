/********************** cmplx.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* make a complex number from two real numbers */
#include "../include/config.h"
#include "../include/complex.h"

complex cmplx( float x, float y )  {
    complex c;
    c.real = x; c.imag = y;
    return(c);
}
