/********************  addmat.c (in su3.a)  *****************************
*									*
*  Add two SU3 matrices 						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void add_su3_matrix(su3_matrix *a,su3_matrix *b,su3_matrix *c)
{
  register int i;
  register int j;
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) {{
        (c -> e)[i][j].real = ((a -> e)[i][j].real + (b -> e)[i][j].real);
        (c -> e)[i][j].imag = ((a -> e)[i][j].imag + (b -> e)[i][j].imag);
      };
    }
}
