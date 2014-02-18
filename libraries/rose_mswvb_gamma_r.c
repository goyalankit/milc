/************* mswvb_gamma_r.c **************************/
/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_swv_by_gamma_r(spin_wilson_vector *src,spin_wilson_vector *dest,int dir)
{
/*color*/
  register int i;
/* row  spin indices*/
  register int s1;
  switch(dir){
    case 0:
{
      for (i = 0; i < 3; i++) 
        for (s1 = 0; s1 < 4; s1++) {{
            (dest -> d)[s1].d[0].c[i].real = (src -> d)[s1].d[3].c[i].imag;
            (dest -> d)[s1].d[0].c[i].imag = -(src -> d)[s1].d[3].c[i].real;
          };
{
            (dest -> d)[s1].d[1].c[i].real = (src -> d)[s1].d[2].c[i].imag;
            (dest -> d)[s1].d[1].c[i].imag = -(src -> d)[s1].d[2].c[i].real;
          };
{
            (dest -> d)[s1].d[2].c[i].real = -(src -> d)[s1].d[1].c[i].imag;
            (dest -> d)[s1].d[2].c[i].imag = (src -> d)[s1].d[1].c[i].real;
          };
{
            (dest -> d)[s1].d[3].c[i].real = -(src -> d)[s1].d[0].c[i].imag;
            (dest -> d)[s1].d[3].c[i].imag = (src -> d)[s1].d[0].c[i].real;
          };
        }
      break; 
    }
    case 1:
{
      for (i = 0; i < 3; i++) 
        for (s1 = 0; s1 < 4; s1++) {{
            (dest -> d)[s1].d[0].c[i].real = -(src -> d)[s1].d[3].c[i].real;
            (dest -> d)[s1].d[0].c[i].imag = -(src -> d)[s1].d[3].c[i].imag;
          };
{
            (dest -> d)[s1].d[1].c[i].real = (src -> d)[s1].d[2].c[i].real;
            (dest -> d)[s1].d[1].c[i].imag = (src -> d)[s1].d[2].c[i].imag;
          };
{
            (dest -> d)[s1].d[2].c[i].real = (src -> d)[s1].d[1].c[i].real;
            (dest -> d)[s1].d[2].c[i].imag = (src -> d)[s1].d[1].c[i].imag;
          };
{
            (dest -> d)[s1].d[3].c[i].real = -(src -> d)[s1].d[0].c[i].real;
            (dest -> d)[s1].d[3].c[i].imag = -(src -> d)[s1].d[0].c[i].imag;
          };
        }
      break; 
    }
    case 2:
{
      for (i = 0; i < 3; i++) 
        for (s1 = 0; s1 < 4; s1++) {{
            (dest -> d)[s1].d[0].c[i].real = (src -> d)[s1].d[2].c[i].imag;
            (dest -> d)[s1].d[0].c[i].imag = -(src -> d)[s1].d[2].c[i].real;
          };
{
            (dest -> d)[s1].d[1].c[i].real = -(src -> d)[s1].d[3].c[i].imag;
            (dest -> d)[s1].d[1].c[i].imag = (src -> d)[s1].d[3].c[i].real;
          };
{
            (dest -> d)[s1].d[2].c[i].real = -(src -> d)[s1].d[0].c[i].imag;
            (dest -> d)[s1].d[2].c[i].imag = (src -> d)[s1].d[0].c[i].real;
          };
{
            (dest -> d)[s1].d[3].c[i].real = (src -> d)[s1].d[1].c[i].imag;
            (dest -> d)[s1].d[3].c[i].imag = -(src -> d)[s1].d[1].c[i].real;
          };
        }
      break; 
    }
    case 3:
{
      for (i = 0; i < 3; i++) 
        for (s1 = 0; s1 < 4; s1++) {{
            (dest -> d)[s1].d[0].c[i].real = (src -> d)[s1].d[2].c[i].real;
            (dest -> d)[s1].d[0].c[i].imag = (src -> d)[s1].d[2].c[i].imag;
          };
{
            (dest -> d)[s1].d[1].c[i].real = (src -> d)[s1].d[3].c[i].real;
            (dest -> d)[s1].d[1].c[i].imag = (src -> d)[s1].d[3].c[i].imag;
          };
{
            (dest -> d)[s1].d[2].c[i].real = (src -> d)[s1].d[0].c[i].real;
            (dest -> d)[s1].d[2].c[i].imag = (src -> d)[s1].d[0].c[i].imag;
          };
{
            (dest -> d)[s1].d[3].c[i].real = (src -> d)[s1].d[1].c[i].real;
            (dest -> d)[s1].d[3].c[i].imag = (src -> d)[s1].d[1].c[i].imag;
          };
        }
      break; 
    }
    case -1:
{
      for (i = 0; i < 3; i++) 
        for (s1 = 0; s1 < 4; s1++) {{
            (dest -> d)[s1].d[0].c[i].real = (src -> d)[s1].d[0].c[i].real;
            (dest -> d)[s1].d[0].c[i].imag = (src -> d)[s1].d[0].c[i].imag;
          };
{
            (dest -> d)[s1].d[1].c[i].real = (src -> d)[s1].d[1].c[i].real;
            (dest -> d)[s1].d[1].c[i].imag = (src -> d)[s1].d[1].c[i].imag;
          };
{
            (dest -> d)[s1].d[2].c[i].real = -(src -> d)[s1].d[2].c[i].real;
            (dest -> d)[s1].d[2].c[i].imag = -(src -> d)[s1].d[2].c[i].imag;
          };
{
            (dest -> d)[s1].d[3].c[i].real = -(src -> d)[s1].d[3].c[i].real;
            (dest -> d)[s1].d[3].c[i].imag = -(src -> d)[s1].d[3].c[i].imag;
          };
        }
      break; 
    }
    default:
{
      printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
    }
  }
}
