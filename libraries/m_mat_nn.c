/*******************  m_mat_nn.c  (in su3.a) ****************************
 *									*
 * void mult_su3_nn( su3_matrix *a,*b,*c )				*
 * matrix multiply, no adjoints 						*
 * C  <-  A*B								*
 */
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

#ifndef FAST
#ifndef NDS

void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    register int i,j,k;
    register complex x,y;

  // Warm up the cache to ensure fair comparison.
  float sum = 0.0f;
  for(i=0;i<3;i++){
      for(j=0;j<3;j++){
          sum += a->e[i][j].real + a->e[i][j].imag;
          sum += b->e[i][j].real + b->e[i][j].imag;
          sum += c->e[i][j].real + c->e[i][j].imag;
      }
  }

    double tick = clock();
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++){
            x.real=x.imag=0.0;
            for(k=0;k<3;k++){
                CMUL( a->e[i][k] , b->e[k][j] , y );
                CSUM( x , y );
            }
            c->e[i][j] = x;
        }
    }
    mat_time += (clock() - tick);
}

#else

void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    register int i,j,k;
    register complex x,y;

    su3_matrix_mod am = su3_matrix_to_su3_matrix_mod(a);
    su3_matrix_mod bm = su3_matrix_to_su3_matrix_mod(b);
    su3_matrix_mod cm = su3_matrix_to_su3_matrix_mod(c);

    double tick = clock();
#pragma simd
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++){
            x.real=x.imag=0.0;
            for(k=0;k<3;k++){
                MCMULJ_( am.real[i][k], am.imag[i][k] , bm.real[k][j], bm.imag[k][j], y );
                CSUM( x , y );
            }

            cm.real[i][j] = x.real;
            cm.imag[i][j] = x.imag;
        }
    }
    mat_time += (clock() - tick);

    *a = su3_matrix_mod_to_su3_matrix(am);
    *b = su3_matrix_mod_to_su3_matrix(bm);
    *c = su3_matrix_mod_to_su3_matrix(cm);
}




#endif

/* "Hand coded" routines, clearer coding is up above */
#else
#ifdef NATIVEDOUBLE   /* RS6000 version */

void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int j;
    register double a0r,a0i,a1r,a1i,a2r,a2i;
    register double b0r,b0i,b1r,b1i,b2r,b2i;

  // Warm up the cache to ensure fair comparison.
  float sum = 0.0f;
  for(i=0;i<3;i++){
      for(j=0;j<3;j++){
          sum += a->e[i][j].real + a->e[i][j].imag;
          sum += b->e[i][j].real + b->e[i][j].imag;
          sum += c->e[i][j].real + c->e[i][j].imag;
      }
  }

    double tick = clock();
    for(j=0;j<3;j++){

        a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
        b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
        a1r=a->e[0][1].real; a1i=a->e[0][1].imag;
        b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
        a2r=a->e[0][2].real; a2i=a->e[0][2].imag;
        b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

        c->e[0][j].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
        c->e[0][j].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

        a0r=a->e[1][0].real; a0i=a->e[1][0].imag;
        b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
        a1r=a->e[1][1].real; a1i=a->e[1][1].imag;
        b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
        a2r=a->e[1][2].real; a2i=a->e[1][2].imag;
        b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

        c->e[1][j].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
        c->e[1][j].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

        a0r=a->e[2][0].real; a0i=a->e[2][0].imag;
        b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
        a1r=a->e[2][1].real; a1i=a->e[2][1].imag;
        b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
        a2r=a->e[2][2].real; a2i=a->e[2][2].imag;
        b2r=b->e[2][j].real; b2i=b->e[2][j].imag;

        c->e[2][j].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
        c->e[2][j].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

    }
    mat_time += (clock() - tick);
}
#else

void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int i,j;
    register float t,ar,ai,br,bi,cr,ci;

  // Warm up the cache to ensure fair comparison.
  float sum = 0.0f;
  for(i=0;i<3;i++){
      for(j=0;j<3;j++){
          sum += a->e[i][j].real + a->e[i][j].imag;
          sum += b->e[i][j].real + b->e[i][j].imag;
          sum += c->e[i][j].real + c->e[i][j].imag;
      }
  }
    double tick = clock();
    for(i=0;i<3;i++)for(j=0;j<3;j++){

        ar=a->e[i][0].real; ai=a->e[i][0].imag;
        br=b->e[0][j].real; bi=b->e[0][j].imag;
        cr=ar*br; t=ai*bi; cr -= t;
        ci=ar*bi; t=ai*br; ci += t;

        ar=a->e[i][1].real; ai=a->e[i][1].imag;
        br=b->e[1][j].real; bi=b->e[1][j].imag;
        t=ar*br; cr += t; t=ai*bi; cr -= t;
        t=ar*bi; ci += t; t=ai*br; ci += t;

        ar=a->e[i][2].real; ai=a->e[i][2].imag;
        br=b->e[2][j].real; bi=b->e[2][j].imag;
        t=ar*br; cr += t; t=ai*bi; cr -= t;
        t=ar*bi; ci += t; t=ai*br; ci += t;

        c->e[i][j].real=cr;
        c->e[i][j].imag=ci;
    }
    mat_time += (clock() - tick);
}
#endif	/* End of "#ifdef NATIVEDOUBLE" */
#endif	/* End of "#ifdef FAST" */
