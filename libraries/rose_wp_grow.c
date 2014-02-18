/***************** wp_grow.c  (in su3.a) **************************/
/* 
  Expand the "Wilson projection" of a Wilson fermion vector.
  (1 +- gamma_j) is a projection operator, and we are given a
  half_wilson_vector which contains the two components of a Wilson
  vector projected out.  This routine reexpands it to a four component
  object.
  usage:  wp_grow(  half_wilson_vector *src, wilson_vector *dest,
        int dir, int sign );
	If dir is one of XUP,YUP,ZUP or TUP, the projection is
	along the eigenvectors with eigenvalue +1, which survive
	multiplcation by (1+gamma[dir]).
	If dir is one of XDOWN,YDOWN,ZDOWN or TDOWN, the projection is
	along the eigenvectors with eigenvalue -1, which survive
	multiplication by (1-gamma[OPP_DIR(dir)]).
	If sign=MINUS reverse the roles of +1 and -1 - in other words
	use -gamma_dir instead of gamma_dir
  Here my eigenvectors are normalized to 2, so for XYZT directions
  I won't explicitely multiply by 2.  In other words, the matrix of
  eigenvectors is sqrt(2) times a unitary matrix, and in reexpanding
  the vector I will multiply by the adjoint of this matrix.
  For UP directions, hvec.h[0] and hvec.h[2] contain the projections
  along the first and second eigenvectors respectively.
  For DOWN directions, hvec.h[0] and hvec.h[2] contain the projections
  along the third and fourth eigenvectors respectively. This results
  in down directions differing from up directions only in the sign of
  the addition.
  Note: wp_shrink( +-dir) followed by wp_grow( +-dir) amounts to multiplication
   by 1+-gamma_dir
 gamma(XUP) 			eigenvectors	eigenvalue
 	    0  0  0  i		( 1, 0, 0,-i)	+1
            0  0  i  0		( 0, 1,-i, 0)	+1
            0 -i  0  0		( 0, 1, 0,+i)	-1
           -i  0  0  0		( 1, 0,+i, 0)	-1
 gamma(YUP)			eigenvectors	eigenvalue
 	    0  0  0 -1		( 1, 0, 0,-1)	+1
            0  0  1  0		( 0, 1, 1, 0)	+1
            0  1  0  0		( 1, 0, 0, 1)	-1
           -1  0  0  0		( 0, 1,-1, 0)	-1
 gamma(ZUP)			eigenvectors	eigenvalue
 	    0  0  i  0		( 1, 0,-i, 0)	+1
            0  0  0 -i		( 0, 1, 0,+i)	+1
           -i  0  0  0		( 1, 0,+i, 0)	-1
            0  i  0  0		( 0, 1, 0,-i)	-1
 gamma(TUP)			eigenvectors	eigenvalue
 	    0  0  1  0		( 1, 0, 1, 0)	+1
            0  0  0  1		( 0, 1, 0, 1)	+1
            1  0  0  0		( 1, 0,-1, 0)	-1
            0  1  0  0		( 0, 1, 0,-1)	-1
 gamma(FIVE) 			eigenvectors	eigenvalue
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void wp_grow(half_wilson_vector *src,wilson_vector *dest,int dir,int sign)
{
/*color*/
  register int i;
/* two ways to get -gamma_dir ! */
  if (sign == -1) 
    dir = (7 - dir);
  switch(dir){
    case 0:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[3].c[i].real = (src -> h)[0].c[i].imag;
          (dest -> d)[3].c[i].imag = -(src -> h)[0].c[i].real;
        };
{
          (dest -> d)[2].c[i].real = (src -> h)[1].c[i].imag;
          (dest -> d)[2].c[i].imag = -(src -> h)[1].c[i].real;
        };
      }
      break; 
    }
    case 7:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[3].c[i].real = -(src -> h)[0].c[i].imag;
          (dest -> d)[3].c[i].imag = (src -> h)[0].c[i].real;
        };
{
          (dest -> d)[2].c[i].real = -(src -> h)[1].c[i].imag;
          (dest -> d)[2].c[i].imag = (src -> h)[1].c[i].real;
        };
      }
      break; 
    }
    case 1:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[3].c[i].real = -(src -> h)[0].c[i].real;
          (dest -> d)[3].c[i].imag = -(src -> h)[0].c[i].imag;
        };
{
          (dest -> d)[2].c[i].real = (src -> h)[1].c[i].real;
          (dest -> d)[2].c[i].imag = (src -> h)[1].c[i].imag;
        };
      }
      break; 
    }
    case 6:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[3].c[i].real = (src -> h)[0].c[i].real;
          (dest -> d)[3].c[i].imag = (src -> h)[0].c[i].imag;
        };
{
          (dest -> d)[2].c[i].real = -(src -> h)[1].c[i].real;
          (dest -> d)[2].c[i].imag = -(src -> h)[1].c[i].imag;
        };
      }
      break; 
    }
    case 2:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[2].c[i].real = (src -> h)[0].c[i].imag;
          (dest -> d)[2].c[i].imag = -(src -> h)[0].c[i].real;
        };
{
          (dest -> d)[3].c[i].real = -(src -> h)[1].c[i].imag;
          (dest -> d)[3].c[i].imag = (src -> h)[1].c[i].real;
        };
      }
      break; 
    }
    case 5:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[2].c[i].real = -(src -> h)[0].c[i].imag;
          (dest -> d)[2].c[i].imag = (src -> h)[0].c[i].real;
        };
{
          (dest -> d)[3].c[i].real = (src -> h)[1].c[i].imag;
          (dest -> d)[3].c[i].imag = -(src -> h)[1].c[i].real;
        };
      }
      break; 
    }
    case 3:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
        (dest -> d)[2].c[i] = (src -> h)[0].c[i];
        (dest -> d)[3].c[i] = (src -> h)[1].c[i];
      }
      break; 
    }
    case 4:
{
      for (i = 0; i < 3; i++) {
        (dest -> d)[0].c[i] = (src -> h)[0].c[i];
        (dest -> d)[1].c[i] = (src -> h)[1].c[i];
{
          (dest -> d)[2].c[i].real = -(src -> h)[0].c[i].real;
          (dest -> d)[2].c[i].imag = -(src -> h)[0].c[i].imag;
        };
{
          (dest -> d)[3].c[i].real = -(src -> h)[1].c[i].real;
          (dest -> d)[3].c[i].imag = -(src -> h)[1].c[i].imag;
        };
      }
      break; 
    }
    default:
{
      printf("BAD CALL TO WP_GROW()\n");
    }
  }
}