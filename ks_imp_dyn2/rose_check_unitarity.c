/*********************** check_unitarity.c ***************************/
/* MIMD version 6 */
/* Claude Bernard, original version */
/* Modified 7/96 DT to quit after finding violation */
/* 9/4/96 DT added Urs' row orthogonality checks, if STRONG defined */
/* 11/2/98 CD fixed return value of max deviation */
/* 01/20/00 UMH combined with Schroedinger functional version */
#define STRONG	/* check row orthogonality as well as norms */
/* Check unitarity of link matrices, quit if not unitary */
#include "generic_includes.h"
#define TOLERANCE (0.0001)
/*#define UNIDEBUG */
float check_su3(su3_matrix *c);

float check_unitarity()
{
  register int i;
  register int dir;
  int ii;
  int jj;
  register site *s;
  register su3_matrix *mat;
  float deviation;
  float max_deviation;
  double av_deviation;
  union __unnamed_class___F0_L26_C2_L70R_variable_declaration__variable_type_f_variable_name_L70R__scope__fval__DELIMITER__L70R_variable_declaration__variable_type_i_variable_name_L70R__scope__ival {
  float fval;
  int ival;}ifval;
  max_deviation = (av_deviation = 0);
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
#ifdef SCHROED_FUN
#else
    for (dir = 0; dir <= 3; dir++) {
#endif
      mat = ((s -> link) + dir);
      deviation = check_su3(mat);
      if (deviation > 0.0001) {
        printf("Unitarity problem on node %d, site %d, dir %d, deviation=%f\n",mynode(),i,dir,deviation);
        printf("SU3 matrix:\n");
        for (ii = 0; ii <= 2; ii++) {
          for (jj = 0; jj <= 2; jj++) {
            printf("%f ",(mat -> e)[ii][jj].real);
            printf("%f ",(mat -> e)[ii][jj].imag);
          }
          printf("\n");
        }
        printf("repeat in hex:\n");
        for (ii = 0; ii <= 2; ii++) {
          for (jj = 0; jj <= 2; jj++) {
            ifval.fval = (mat -> e)[ii][jj].real;
            printf("%08x ",ifval.ival);
            ifval.fval = (mat -> e)[ii][jj].imag;
            printf("%08x ",ifval.ival);
          }
          printf("\n");
        }
        printf("  \n \n");
        fflush(stdout);
        terminate(1);
      }
      if (max_deviation < deviation) 
        max_deviation = deviation;
      av_deviation += (deviation * deviation);
    }
  }
  av_deviation = sqrt((av_deviation / (4 * i)));
#ifdef UNIDEBUG
#endif
  if (max_deviation > 0.0001) 
    printf("Unitarity problem on node %d, maximum deviation=%f\n",mynode(),max_deviation);
  return max_deviation;
/*check_unitarity() */
}

float check_su3(su3_matrix *c)
{
  register float ar;
  register float ai;
  register float ari;
  register float max;
  register int i;
/* first normalize row */
  for (((i = 0) , (max = 0.0)); i < 3; ++i) {
/* sum of squares of row */
    ar = (((((((c -> e)[i][0].real * (c -> e)[i][0].real) + ((c -> e)[i][0].imag * (c -> e)[i][0].imag)) + ((c -> e)[i][1].real * (c -> e)[i][1].real)) + ((c -> e)[i][1].imag * (c -> e)[i][1].imag)) + ((c -> e)[i][2].real * (c -> e)[i][2].real)) + ((c -> e)[i][2].imag * (c -> e)[i][2].imag));
    ar = (fabs((sqrt(((double )ar)) - 1.)));
    if (max < ar) 
      max = ar;
  }
#ifdef STRONG
/* Test orthogonality of row 0 and row 1 */
/* real part of 0 dot 1 */
  ar = (((((((c -> e)[0][0].real * (c -> e)[1][0].real) + ((c -> e)[0][0].imag * (c -> e)[1][0].imag)) + ((c -> e)[0][1].real * (c -> e)[1][1].real)) + ((c -> e)[0][1].imag * (c -> e)[1][1].imag)) + ((c -> e)[0][2].real * (c -> e)[1][2].real)) + ((c -> e)[0][2].imag * (c -> e)[1][2].imag));
/* imag part of 0 dot 1 */
  ai = (((((((c -> e)[0][0].real * (c -> e)[1][0].imag) - ((c -> e)[0][0].imag * (c -> e)[1][0].real)) + ((c -> e)[0][1].real * (c -> e)[1][1].imag)) - ((c -> e)[0][1].imag * (c -> e)[1][1].real)) + ((c -> e)[0][2].real * (c -> e)[1][2].imag)) - ((c -> e)[0][2].imag * (c -> e)[1][2].real));
  ari = (sqrt(((double )((ar * ar) + (ai * ai)))));
  if (max < ari) 
    max = ari;
/* Test orthogonality of row 0 and row 2 */
/* real part of 0 dot 1 */
  ar = (((((((c -> e)[0][0].real * (c -> e)[2][0].real) + ((c -> e)[0][0].imag * (c -> e)[2][0].imag)) + ((c -> e)[0][1].real * (c -> e)[2][1].real)) + ((c -> e)[0][1].imag * (c -> e)[2][1].imag)) + ((c -> e)[0][2].real * (c -> e)[2][2].real)) + ((c -> e)[0][2].imag * (c -> e)[2][2].imag));
/* imag part of 0 dot 1 */
  ai = (((((((c -> e)[0][0].real * (c -> e)[2][0].imag) - ((c -> e)[0][0].imag * (c -> e)[2][0].real)) + ((c -> e)[0][1].real * (c -> e)[2][1].imag)) - ((c -> e)[0][1].imag * (c -> e)[2][1].real)) + ((c -> e)[0][2].real * (c -> e)[2][2].imag)) - ((c -> e)[0][2].imag * (c -> e)[2][2].real));
  ari = (sqrt(((double )((ar * ar) + (ai * ai)))));
  if (max < ari) 
    max = ari;
/* Test orthogonality of row 1 and row 2 */
/* real part of 0 dot 1 */
  ar = (((((((c -> e)[1][0].real * (c -> e)[2][0].real) + ((c -> e)[1][0].imag * (c -> e)[2][0].imag)) + ((c -> e)[1][1].real * (c -> e)[2][1].real)) + ((c -> e)[1][1].imag * (c -> e)[2][1].imag)) + ((c -> e)[1][2].real * (c -> e)[2][2].real)) + ((c -> e)[1][2].imag * (c -> e)[2][2].imag));
/* imag part of 0 dot 1 */
  ai = (((((((c -> e)[1][0].real * (c -> e)[2][0].imag) - ((c -> e)[1][0].imag * (c -> e)[2][0].real)) + ((c -> e)[1][1].real * (c -> e)[2][1].imag)) - ((c -> e)[1][1].imag * (c -> e)[2][1].real)) + ((c -> e)[1][2].real * (c -> e)[2][2].imag)) - ((c -> e)[1][2].imag * (c -> e)[2][2].real));
  ari = (sqrt(((double )((ar * ar) + (ai * ai)))));
  if (max < ari) 
    max = ari;
#endif /*STRONG*/
  return max;
/* check_su3 */
}
