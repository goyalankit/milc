/*********************** reunitarize2.c ***************************/
/* MIMD version 6 */
/* reunitarize the link matrices */
/* Modifications ...
   09/03/96 Incorporates unitarity checking C.D.
   01/20/00 combined with Schroedinger functional version - UMH
   03/28/00 Modified reunit_su3 to allow calling program to handle errors.
*/
#include "generic_includes.h"
#define TOLERANCE (0.0001)
#define MAXERRCOUNT 100
/**#define UNIDEBUG**/
float max_deviation;
double av_deviation;
/* canopy qcdlib code - stolen, of course */
#define fixsu3(matrix) \
 { \
    bj0r = (*matrix).e[0][0].real; \
    bj0i = (*matrix).e[0][0].imag; \
    bj1r = (*matrix).e[0][1].real; \
    bj1i = (*matrix).e[0][1].imag; \
    bj2r = (*matrix).e[0][2].real; \
    bj2i = (*matrix).e[0][2].imag; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = bj1r*ar - bj1i*ai; \
    ti = bj1r*ai + bj1i*ar; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = tr - bj2r*ar + bj2i*ai; \
    ti = ti - bj2r*ai - bj2i*ar; \
    (*matrix).e[2][0].real = tr; \
    (*matrix).e[2][0].imag = -ti; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = bj2r*ar - bj2i*ai; \
    ti = bj2r*ai + bj2i*ar; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = tr - bj0r*ar + bj0i*ai; \
    ti = ti - bj0r*ai - bj0i*ar; \
    (*matrix).e[2][1].real = tr; \
    (*matrix).e[2][1].imag = -ti; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = bj0r*ar - bj0i*ai; \
    ti = bj0r*ai + bj0i*ar; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = tr - bj1r*ar + bj1i*ai; \
    ti = ti - bj1r*ai - bj1i*ar; \
    (*matrix).e[2][2].real = tr; \
    (*matrix).e[2][2].imag = -ti; \
 } /* define fixsu3 */

int check_deviation(float deviation)
{
  if (max_deviation < deviation) 
    max_deviation = deviation;
  av_deviation += (deviation * deviation);
  if (deviation > 0.0001) {
    return 1;
  }
  else 
    return 0;
/* check_deviation */
}

void reunit_report_problem_matrix(su3_matrix *mat,int i,int dir)
{
  int ii;
  int jj;
  union __unnamed_class___F0_L78_C3_L69R__L70R__scope____SgSS2___variable_declaration__variable_type_f_variable_name_L69R__L70R__scope____SgSS2____scope__fval__DELIMITER__L69R__L70R__scope____SgSS2___variable_declaration__variable_type_i_variable_name_L69R__L70R__scope____SgSS2____scope__ival {
  float fval;
  int ival;}ifval;
  printf("Unitarity problem on node %d, site %d, dir %d tolerance=%e\n",mynode(),i,dir,0.0001);
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
/* reunit_report_problem_matrix */
}

int reunit_su3(su3_matrix *c)
{
  register float bj0r;
  register float bj0i;
  register float bj1r;
  register float bj1i;
  register float bj2r;
  register float bj2i;
  register float c0r;
  register float c0i;
  register float c1r;
  register float c1i;
  register float c2r;
  register float c2i;
  register float ar;
  register float ai;
  register float tr;
  register float ti;
  float deviation;
  int errors;
  errors = 0;
/* first normalize row 0 */
/* sum of squares of row */
  ar = (((((((c -> e)[0][0].real * (c -> e)[0][0].real) + ((c -> e)[0][0].imag * (c -> e)[0][0].imag)) + ((c -> e)[0][1].real * (c -> e)[0][1].real)) + ((c -> e)[0][1].imag * (c -> e)[0][1].imag)) + ((c -> e)[0][2].real * (c -> e)[0][2].real)) + ((c -> e)[0][2].imag * (c -> e)[0][2].imag));
  deviation = (fabs((ar - 1.)));
  errors += check_deviation(deviation);
/* used to normalize row */
  ar = (1. / sqrt(((double )ar)));
  (c -> e)[0][0].real *= ar;
  (c -> e)[0][0].imag *= ar;
  (c -> e)[0][1].real *= ar;
  (c -> e)[0][1].imag *= ar;
  (c -> e)[0][2].real *= ar;
  (c -> e)[0][2].imag *= ar;
/* now make row 1 orthogonal to row 0 */
/* real part of 0 dot 1 */
  ar = (((((((c -> e)[0][0].real * (c -> e)[1][0].real) + ((c -> e)[0][0].imag * (c -> e)[1][0].imag)) + ((c -> e)[0][1].real * (c -> e)[1][1].real)) + ((c -> e)[0][1].imag * (c -> e)[1][1].imag)) + ((c -> e)[0][2].real * (c -> e)[1][2].real)) + ((c -> e)[0][2].imag * (c -> e)[1][2].imag));
/* imag part of 0 dot 1 */
  ai = (((((((c -> e)[0][0].real * (c -> e)[1][0].imag) - ((c -> e)[0][0].imag * (c -> e)[1][0].real)) + ((c -> e)[0][1].real * (c -> e)[1][1].imag)) - ((c -> e)[0][1].imag * (c -> e)[1][1].real)) + ((c -> e)[0][2].real * (c -> e)[1][2].imag)) - ((c -> e)[0][2].imag * (c -> e)[1][2].real));
  deviation = ((ar * ar) + (ai * ai));
  errors += check_deviation(deviation);
/* row 1 -= a * row 0 */
  (c -> e)[1][0].real -= ((ar * (c -> e)[0][0].real) - (ai * (c -> e)[0][0].imag));
  (c -> e)[1][0].imag -= ((ar * (c -> e)[0][0].imag) + (ai * (c -> e)[0][0].real));
  (c -> e)[1][1].real -= ((ar * (c -> e)[0][1].real) - (ai * (c -> e)[0][1].imag));
  (c -> e)[1][1].imag -= ((ar * (c -> e)[0][1].imag) + (ai * (c -> e)[0][1].real));
  (c -> e)[1][2].real -= ((ar * (c -> e)[0][2].real) - (ai * (c -> e)[0][2].imag));
  (c -> e)[1][2].imag -= ((ar * (c -> e)[0][2].imag) + (ai * (c -> e)[0][2].real));
/* now normalize row 1 */
/* sum of squares of row */
  ar = (((((((c -> e)[1][0].real * (c -> e)[1][0].real) + ((c -> e)[1][0].imag * (c -> e)[1][0].imag)) + ((c -> e)[1][1].real * (c -> e)[1][1].real)) + ((c -> e)[1][1].imag * (c -> e)[1][1].imag)) + ((c -> e)[1][2].real * (c -> e)[1][2].real)) + ((c -> e)[1][2].imag * (c -> e)[1][2].imag));
  deviation = (fabs((ar - 1.)));
  errors += check_deviation(deviation);
/* used to normalize row */
  ar = (1. / sqrt(((double )ar)));
  (c -> e)[1][0].real *= ar;
  (c -> e)[1][0].imag *= ar;
  (c -> e)[1][1].real *= ar;
  (c -> e)[1][1].imag *= ar;
  (c -> e)[1][2].real *= ar;
  (c -> e)[1][2].imag *= ar;
/* Save for checking */
  c0r = (c -> e)[2][0].real;
  c0i = (c -> e)[2][0].imag;
  c1r = (c -> e)[2][1].real;
  c1i = (c -> e)[2][1].imag;
  c2r = (c -> e)[2][2].real;
  c2i = (c -> e)[2][2].imag;
/* reconstruct row 2 */
{
    bj0r = (c -> e)[0][0].real;
    bj0i = (c -> e)[0][0].imag;
    bj1r = (c -> e)[0][1].real;
    bj1i = (c -> e)[0][1].imag;
    bj2r = (c -> e)[0][2].real;
    bj2i = (c -> e)[0][2].imag;
    ar = (c -> e)[1][2].real;
    ai = (c -> e)[1][2].imag;
    tr = ((bj1r * ar) - (bj1i * ai));
    ti = ((bj1r * ai) + (bj1i * ar));
    ar = (c -> e)[1][1].real;
    ai = (c -> e)[1][1].imag;
    tr = ((tr - (bj2r * ar)) + (bj2i * ai));
    ti = ((ti - (bj2r * ai)) - (bj2i * ar));
    (c -> e)[2][0].real = tr;
    (c -> e)[2][0].imag = -ti;
    ar = (c -> e)[1][0].real;
    ai = (c -> e)[1][0].imag;
    tr = ((bj2r * ar) - (bj2i * ai));
    ti = ((bj2r * ai) + (bj2i * ar));
    ar = (c -> e)[1][2].real;
    ai = (c -> e)[1][2].imag;
    tr = ((tr - (bj0r * ar)) + (bj0i * ai));
    ti = ((ti - (bj0r * ai)) - (bj0i * ar));
    (c -> e)[2][1].real = tr;
    (c -> e)[2][1].imag = -ti;
    ar = (c -> e)[1][1].real;
    ai = (c -> e)[1][1].imag;
    tr = ((bj0r * ar) - (bj0i * ai));
    ti = ((bj0r * ai) + (bj0i * ar));
    ar = (c -> e)[1][0].real;
    ai = (c -> e)[1][0].imag;
    tr = ((tr - (bj1r * ar)) + (bj1i * ai));
    ti = ((ti - (bj1r * ai)) - (bj1i * ar));
    (c -> e)[2][2].real = tr;
    (c -> e)[2][2].imag = -ti;
  };
/* Now check deviation */
  ar = (((((((c0r - (c -> e)[2][0].real) * (c0r - (c -> e)[2][0].real)) + ((c0i - (c -> e)[2][0].imag) * (c0i - (c -> e)[2][0].imag))) + ((c1r - (c -> e)[2][1].real) * (c1r - (c -> e)[2][1].real))) + ((c1i - (c -> e)[2][1].imag) * (c1i - (c -> e)[2][1].imag))) + ((c2r - (c -> e)[2][2].real) * (c2r - (c -> e)[2][2].real))) + ((c2i - (c -> e)[2][2].imag) * (c2i - (c -> e)[2][2].imag)));
  deviation = ar;
  errors += check_deviation(deviation);
  return errors;
/* reunit_su3 */
}

void reunitarize()
{
  register su3_matrix *mat;
  register int i;
  register int dir;
  register site *s;
  int errcount = 0;
  int errors;
  max_deviation = 0.0;
  av_deviation = 0.0;
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
#ifdef SCHROED_FUN
#else
    for (dir = 0; dir <= 3; dir++) {
#endif
      mat = ((s -> link) + dir);
      errors = reunit_su3(mat);
      errcount += errors;
      if (errors != 0) 
        reunit_report_problem_matrix(mat,i,dir);
      if (errcount > 100) {
        printf("Unitarity error count exceeded.\n");
        terminate(1);
      }
    }
  }
#ifdef UNIDEBUG
#endif
  if (max_deviation > 0.0001) {
    printf("reunitarize: Node %d unitarity problem, maximum deviation=%e\n",mynode(),max_deviation);
    errcount++;
    if (errcount > 100) {
      printf("Unitarity error count exceeded.\n");
      terminate(1);
    }
  }
/* reunitarize2 */
}
