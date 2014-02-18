/*****************  wp_shrink4.c  (in su3.a) ****************************
*									*
* Shrink a wilson vector in four directions, producing four		*
*  half_wilson_vectors.							*
* void wp_shrink_4dir(  wilson_vector *a,  half_wilson_vector *b1,	*
*       half_wilson_vector *b2, half_wilson_vector *b3,			*
*       half_wilson_vector *b4, int sign );				*
* B1 <- (1 +- gamma_x)A,, projection					*
*  argument "sign" is sign of gamma matrix.				*
*  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"
#ifndef FAST   /* "FAST", or IBM RS6000 version inlines calls */

void wp_shrink_4dir(wilson_vector *a,half_wilson_vector *b1,half_wilson_vector *b2,half_wilson_vector *b3,half_wilson_vector *b4,int sign)
{
  wp_shrink(a,b1,0,sign);
  wp_shrink(a,b2,1,sign);
  wp_shrink(a,b3,2,sign);
  wp_shrink(a,b4,3,sign);
}
#else   /* "FAST" code inlines calls */
/*color*/
/*    wp_shrink( a,b1,XUP,sign); */
/* case XUP: */
/* case XDOWN: */
/*    wp_shrink( a,b2,YUP,sign); */
/* case YUP: */
/* case YDOWN: */
/*    wp_shrink( a,b3,ZUP,sign); */
/* case ZUP: */
/* case ZDOWN: */
/*    wp_shrink( a,b4,TUP,sign); */
/* case TUP: */
/* case TDOWN: */
#endif /* "ifndef FAST */
