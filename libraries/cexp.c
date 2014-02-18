/********************** cexp.c (in complex.a) **********************/
/* MIMD version 6 */
/* Subroutines for operations on complex numbers */
/* complex exponential */
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"

complex cexp( complex *a ){
    complex c;
    float mag;
    mag = (float)exp( (double)(*a).real );
    c.real = mag*(float)cos( (double)(*a).imag );
    c.imag = mag*(float)sin( (double)(*a).imag );
    return(c);
}
