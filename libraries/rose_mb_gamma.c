/************* mb_gamma.c  (in su3.a) **************************/
/* 
  Multiply a Wilson vector by a gamma matrix
  usage:  mult_by_gamma( wilson_vector *src, wilson_vector *dest, int dir )
	dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0
 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0
 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0
 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0
 gamma(FIVE) 
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_by_gamma(wilson_vector *src,wilson_vector *dest,int dir)
{
/*color*/
  register int i;
  switch(dir){
    case 0:
{
      for (i = 0; i < 3; i++) {{
          (dest -> d)[0].c[i].real = -(src -> d)[3].c[i].imag;
          (dest -> d)[0].c[i].imag = (src -> d)[3].c[i].real;
        };
{
          (dest -> d)[1].c[i].real = -(src -> d)[2].c[i].imag;
          (dest -> d)[1].c[i].imag = (src -> d)[2].c[i].real;
        };
{
          (dest -> d)[2].c[i].real = (src -> d)[1].c[i].imag;
          (dest -> d)[2].c[i].imag = -(src -> d)[1].c[i].real;
        };
{
          (dest -> d)[3].c[i].real = (src -> d)[0].c[i].imag;
          (dest -> d)[3].c[i].imag = -(src -> d)[0].c[i].real;
        };
      }
      break; 
    }
    case 1:
{
      for (i = 0; i < 3; i++) {{
          (dest -> d)[0].c[i].real = -(src -> d)[3].c[i].real;
          (dest -> d)[0].c[i].imag = -(src -> d)[3].c[i].imag;
        };
{
          (dest -> d)[1].c[i].real = (src -> d)[2].c[i].real;
          (dest -> d)[1].c[i].imag = (src -> d)[2].c[i].imag;
        };
{
          (dest -> d)[2].c[i].real = (src -> d)[1].c[i].real;
          (dest -> d)[2].c[i].imag = (src -> d)[1].c[i].imag;
        };
{
          (dest -> d)[3].c[i].real = -(src -> d)[0].c[i].real;
          (dest -> d)[3].c[i].imag = -(src -> d)[0].c[i].imag;
        };
      }
      break; 
    }
    case 2:
{
      for (i = 0; i < 3; i++) {{
          (dest -> d)[0].c[i].real = -(src -> d)[2].c[i].imag;
          (dest -> d)[0].c[i].imag = (src -> d)[2].c[i].real;
        };
{
          (dest -> d)[1].c[i].real = (src -> d)[3].c[i].imag;
          (dest -> d)[1].c[i].imag = -(src -> d)[3].c[i].real;
        };
{
          (dest -> d)[2].c[i].real = (src -> d)[0].c[i].imag;
          (dest -> d)[2].c[i].imag = -(src -> d)[0].c[i].real;
        };
{
          (dest -> d)[3].c[i].real = -(src -> d)[1].c[i].imag;
          (dest -> d)[3].c[i].imag = (src -> d)[1].c[i].real;
        };
      }
      break; 
    }
    case 3:
{
      for (i = 0; i < 3; i++) {{
          (dest -> d)[0].c[i].real = (src -> d)[2].c[i].real;
          (dest -> d)[0].c[i].imag = (src -> d)[2].c[i].imag;
        };
{
          (dest -> d)[1].c[i].real = (src -> d)[3].c[i].real;
          (dest -> d)[1].c[i].imag = (src -> d)[3].c[i].imag;
        };
{
          (dest -> d)[2].c[i].real = (src -> d)[0].c[i].real;
          (dest -> d)[2].c[i].imag = (src -> d)[0].c[i].imag;
        };
{
          (dest -> d)[3].c[i].real = (src -> d)[1].c[i].real;
          (dest -> d)[3].c[i].imag = (src -> d)[1].c[i].imag;
        };
      }
      break; 
    }
    case -1:
{
      for (i = 0; i < 3; i++) {{
          (dest -> d)[0].c[i].real = (src -> d)[0].c[i].real;
          (dest -> d)[0].c[i].imag = (src -> d)[0].c[i].imag;
        };
{
          (dest -> d)[1].c[i].real = (src -> d)[1].c[i].real;
          (dest -> d)[1].c[i].imag = (src -> d)[1].c[i].imag;
        };
{
          (dest -> d)[2].c[i].real = -(src -> d)[2].c[i].real;
          (dest -> d)[2].c[i].imag = -(src -> d)[2].c[i].imag;
        };
{
          (dest -> d)[3].c[i].real = -(src -> d)[3].c[i].real;
          (dest -> d)[3].c[i].imag = -(src -> d)[3].c[i].imag;
        };
      }
      break; 
    }
    default:
{
      printf("BAD CALL TO MULT_BY_GAMMA()\n");
    }
  }
}
