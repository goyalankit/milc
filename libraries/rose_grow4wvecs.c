/*****************  grow4wvecs.c  (in su3.a) ****************************
*									*
*  If sum=0,								*
*  Grow and add four wilson_vectors 					*
*  If sum=1,								*
*  Grow and sum four wilson_vectors to another wilson_vector		*
* void grow_four_wvecs(a,b1,b2,b3,b4,sign,sum)				*
* wilson_vector *a; half_wilson_vector *b1,*b2,*b3,*b4;			*
* int sign,sum;								*
* A  <-  B1 + B2 + B3 + B4   or						*
* A  <-  A + B1 + B2 + B3 + B4						*
* B1 is expanded using gamma_x, B2 using gamma_y, etc. 			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"
/* grow and sum four wilson_vectors */
#ifndef FAST

void grow_add_four_wvecs(wilson_vector *a,half_wilson_vector *b1,half_wilson_vector *b2,half_wilson_vector *b3,half_wilson_vector *b4,int sign,int sum)
{
  if (sum == 0) 
    wp_grow(b1,a,0,sign);
  else 
    wp_grow_add(b1,a,0,sign);
  wp_grow_add(b2,a,1,sign);
  wp_grow_add(b3,a,2,sign);
  wp_grow_add(b4,a,3,sign);
}
#else  /* "FAST" code has wp_grow_add inlined */
/* For the RS6000 */
/* wp_grow( b1,a,XUP,sign); */
/* case XUP: */
/* case XDOWN: */
/*wp_grow_add( b1,a,XUP,sign); */
/* case XUP: */
/* case XDOWN: */
/* wp_grow_add( b2,a,YUP,sign); */
/* case YUP: */
/* case YDOWN: */
/* wp_grow_add( b3,a,ZUP,sign); */
/* case ZUP: */
/* case ZDOWN:*/
/* wp_grow_add( b4,a,TUP,sign); */
/* case TUP: */
/* case TDOWN: */
#endif /* "#ifndef FAST"*/
