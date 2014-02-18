#ifndef MATRIX_H
#define MATRIX_H

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

static inline complex cmplx( float x, float y )  {
    complex c;
    c.real = x; c.imag = y;
    return(c);
}

static inline void add_su3_matrix( su3_matrix *a, su3_matrix *b, su3_matrix *c ) {
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CADD( a->e[i][j], b->e[i][j], c->e[i][j] );
    }
}

static inline void scalar_mult_add_su3_matrix(su3_matrix *a,su3_matrix *b,float s, su3_matrix *c){

register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
	c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
    }
}

static inline void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
  int j;
  register float a0r,a0i,a1r,a1i,a2r,a2i;
  register float b0r,b0i,b1r,b1i,b2r,b2i;
  
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
}

static inline void mult_adj_su3_mat_hwvec( su3_matrix *mat,
       half_wilson_vector *src, half_wilson_vector *dest ){

  register float a0r,a0i,a1r,a1i,a2r,a2i;
  register float b0r,b0i,b1r,b1i,b2r,b2i;

/*    mult_adj_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */
  
  a0r=mat->e[0][0].real;   a0i=mat->e[0][0].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][0].real;   a1i=mat->e[1][0].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][0].real;   a2i=mat->e[2][0].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][1].real;   a0i=mat->e[0][1].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;  
  a1r=mat->e[1][1].real;   a1i=mat->e[1][1].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][1].real;   a2i=mat->e[2][1].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][2].real;   a0i=mat->e[0][2].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][2].real;   a1i=mat->e[1][2].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][2].real;   a2i=mat->e[2][2].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;


/*    mult_adj_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */

  a0r=mat->e[0][0].real;   a0i=mat->e[0][0].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][0].real;   a1i=mat->e[1][0].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][0].real;   a2i=mat->e[2][0].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][1].real;   a0i=mat->e[0][1].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;  
  a1r=mat->e[1][1].real;   a1i=mat->e[1][1].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][1].real;   a2i=mat->e[2][1].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][2].real;   a0i=mat->e[0][2].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][2].real;   a1i=mat->e[1][2].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][2].real;   a2i=mat->e[2][2].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;

}

#if 0
static inline void mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
int i,j;
register float t,ar,ai,br,bi,cr,ci;
    for(i=0;i<3;i++)for(j=0;j<3;j++){

	ar=a->e[i][0].real; ai=a->e[i][0].imag;
	br=b->e[j][0].real; bi=b->e[j][0].imag;
	cr=ar*br; t=ai*bi; cr += t;
	ci=ai*br; t=ar*bi; ci -= t;

	ar=a->e[i][1].real; ai=a->e[i][1].imag;
	br=b->e[j][1].real; bi=b->e[j][1].imag;
	t=ar*br; cr += t; t=ai*bi; cr += t;
	t=ar*bi; ci -= t; t=ai*br; ci += t;

	ar=a->e[i][2].real; ai=a->e[i][2].imag;
	br=b->e[j][2].real; bi=b->e[j][2].imag;
	t=ar*br; cr += t; t=ai*bi; cr += t;
	t=ar*bi; ci -= t; t=ai*br; ci += t;

	c->e[i][j].real=cr;
	c->e[i][j].imag=ci;
    }
}
#endif


static inline void su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c ){

  register int i,j;
  register float  ar,ai,br,bi;

    for(i=0;i<3;i++){
	ar=a->c[i].real;  ai=a->c[i].imag;
	for(j=0;j<3;j++){
	    br=b->c[j].real;  bi=b->c[j].imag;
	    c->e[i][j].real = ar*br + ai*bi;
	    c->e[i][j].imag = ai*br - ar*bi;
	}
    }
}

//#ifndef NDS
/*
 * NDS: instrumented code START
 *
 *
static inline void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
  int j;

  register float a0r,a0i,a1r,a1i,a2r,a2i;
  register float b0r,b0i,b1r,b1i,b2r,b2i;

  for(j=0;j<3;j++){
    
    a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][0].real; a1i=a->e[1][0].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][0].real; a2i=a->e[2][0].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[0][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[0][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    a0r=a->e[0][1].real; a0i=a->e[0][1].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][1].real; a1i=a->e[1][1].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][1].real; a2i=a->e[2][1].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[1][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[1][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    a0r=a->e[0][2].real; a0i=a->e[0][2].imag;
    b0r=b->e[0][j].real; b0i=b->e[0][j].imag;
    a1r=a->e[1][2].real; a1i=a->e[1][2].imag;
    b1r=b->e[1][j].real; b1i=b->e[1][j].imag;
    a2r=a->e[2][2].real; a2i=a->e[2][2].imag;
    b2r=b->e[2][j].real; b2i=b->e[2][j].imag;
    
    c->e[2][j].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
    c->e[2][j].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
    }
}
*/
//#else


//Find:      (.{3})=(\w)->e(\[[\d]\])(\[(.*)\]).(imag)
//Replace:   $1=$2m.imag$3$4

static inline su3_matrix_mod su3_matrix_to_su3_matrix_mod( su3_matrix *a){
    int i, j;
    su3_matrix_mod s3mm;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            s3mm.real[i][j] = a->e[i][j].real;
            s3mm.imag[i][j] = a->e[i][j].imag;
        }
    }
    return s3mm;
}

static inline su3_matrix su3_matrix_mod_to_su3_matrix( su3_matrix_mod s3mm){
    int i, j;
    su3_matrix *a;
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
             a->e[i][j].real = s3mm.real[i][j];
             a->e[i][j].imag = s3mm.imag[i][j];
        }
    }
    return *a;
}

static inline void mult_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
    int j;

    su3_matrix_mod am = su3_matrix_to_su3_matrix_mod(a);
    su3_matrix_mod bm = su3_matrix_to_su3_matrix_mod(b);
    su3_matrix_mod cm = su3_matrix_to_su3_matrix_mod(c);

    register float a0r,a0i,a1r,a1i,a2r,a2i;
    register float b0r,b0i,b1r,b1i,b2r,b2i;
    
    for(j=0;j<3;j++){

        a0r=am.real[0][0]; a0i=am.imag[0][0];
        b0r=bm.real[0][j]; b0i=bm.imag[0][j];
        a1r=am.real[1][0]; a1i=am.imag[1][0];
        b1r=bm.real[1][j]; b1i=bm.imag[1][j];
        a2r=am.real[2][0]; a2i=am.imag[2][0];
        b2r=bm.real[2][j]; b2i=bm.imag[2][j];

        cm.real[0][j] = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
        cm.imag[0][j] = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;

        a0r=am.real[0][1]; a0i=am.imag[0][1];
        b0r=bm.real[0][j]; b0i=bm.imag[0][j];
        a1r=am.real[1][1]; a1i=am.imag[1][1];
        b1r=bm.real[1][j]; b1i=bm.imag[1][j];
        a2r=am.real[2][1]; a2i=am.imag[2][1];
        b2r=bm.real[2][j]; b2i=bm.imag[2][j];


        cm.real[1][j] = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
        cm.imag[1][j] = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;

        a0r=am.real[0][2]; a0i=am.imag[0][2];
        b0r=bm.real[0][j]; b0i=bm.imag[0][j];
        a1r=am.real[1][2]; a1i=am.imag[1][2];
        b1r=bm.real[1][j]; b1i=bm.imag[1][j];
        a2r=am.real[2][2]; a2i=am.imag[2][2];
        b2r=bm.real[2][j]; b2i=bm.imag[2][j];


        cm.real[2][j] = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
        cm.imag[2][j] = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
        
        *a = su3_matrix_mod_to_su3_matrix(am);
        *b = su3_matrix_mod_to_su3_matrix(bm);
        *c = su3_matrix_mod_to_su3_matrix(cm);
    }
}

//#endif /* End of "ifdef NDS" */

/* NDS END */


static inline void mult_su3_mat_hwvec( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest ){

  register float a0r,a0i,a1r,a1i,a2r,a2i;
  register float b0r,b0i,b1r,b1i,b2r,b2i;
  
/*    mult_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

/*    mult_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

}

static inline void scalar_mult_add_su3_vector(su3_vector *a, su3_vector *b, float s,
	su3_vector *c){

  register int i;
  for(i=0;i<3;i++){
    c->c[i].real = a->c[i].real + s*b->c[i].real;
    c->c[i].imag = a->c[i].imag + s*b->c[i].imag;
  }
}

static inline void su3_adjoint( su3_matrix *a, su3_matrix *b ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	CONJG( a->e[j][i], b->e[i][j] );
    }
}

static inline void scalar_mult_su3_matrix( su3_matrix *a, float s, su3_matrix *b ){

  register float ss;

  ss = s;

  b->e[0][0].real = ss*a->e[0][0].real;
  b->e[0][0].imag = ss*a->e[0][0].imag;
  b->e[0][1].real = ss*a->e[0][1].real;
  b->e[0][1].imag = ss*a->e[0][1].imag;
  b->e[0][2].real = ss*a->e[0][2].real;
  b->e[0][2].imag = ss*a->e[0][2].imag;

  b->e[1][0].real = ss*a->e[1][0].real;
  b->e[1][0].imag = ss*a->e[1][0].imag;
  b->e[1][1].real = ss*a->e[1][1].real;
  b->e[1][1].imag = ss*a->e[1][1].imag;
  b->e[1][2].real = ss*a->e[1][2].real;
  b->e[1][2].imag = ss*a->e[1][2].imag;

  b->e[2][0].real = ss*a->e[2][0].real;
  b->e[2][0].imag = ss*a->e[2][0].imag;
  b->e[2][1].real = ss*a->e[2][1].real;
  b->e[2][1].imag = ss*a->e[2][1].imag;
  b->e[2][2].real = ss*a->e[2][2].real;
  b->e[2][2].imag = ss*a->e[2][2].imag;
}

static inline void su3mat_copy( su3_matrix *a, su3_matrix *b ){
register int i,j;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
	b->e[i][j].real = a->e[i][j].real;
	b->e[i][j].imag = a->e[i][j].imag;
    }
}

static inline void mult_su3_mat_vec_sum_4dir(  su3_matrix *a, su3_vector *b0,
       su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c  ){

  register int n;
  register float  c0r,c0i,c1r,c1i,c2r,c2i;
  register float  br,bi,a0,a1,a2;
  register su3_matrix *mat;
  register su3_vector *b;

  c0r = c0i = c1r = c1i = c2r = c2i = 0.0;
  mat = a;

  for(n=0;n<4;n++,mat++){

  switch(n){
    case(0): b=b0; break;
    case(1): b=b1; break;
    case(2): b=b2; break;
    case(3): b=b3; break;
  }

  br=b->c[0].real;    bi=b->c[0].imag;
  a0=mat->e[0][0].real;
  a1=mat->e[1][0].real;
  a2=mat->e[2][0].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=mat->e[0][0].imag;
  a1=mat->e[1][0].imag;
  a2=mat->e[2][0].imag;

  c0r -= a0*bi;
  c1r -= a1*bi;
  c2r -= a2*bi;
  c0i += a0*br;
  c1i += a1*br;
  c2i += a2*br;

  br=b->c[1].real;    bi=b->c[1].imag;
  a0=mat->e[0][1].real;
  a1=mat->e[1][1].real;
  a2=mat->e[2][1].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=mat->e[0][1].imag;
  a1=mat->e[1][1].imag;
  a2=mat->e[2][1].imag;

  c0r -= a0*bi;
  c1r -= a1*bi;
  c2r -= a2*bi;
  c0i += a0*br;
  c1i += a1*br;
  c2i += a2*br;

  br=b->c[2].real;    bi=b->c[2].imag;
  a0=mat->e[0][2].real;
  a1=mat->e[1][2].real;
  a2=mat->e[2][2].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=mat->e[0][2].imag;
  a1=mat->e[1][2].imag;
  a2=mat->e[2][2].imag;

  c0r -= a0*bi;
  c1r -= a1*bi;
  c2r -= a2*bi;
  c0i += a0*br;
  c1i += a1*br;
  c2i += a2*br;

  }

  c->c[0].real = c0r;
  c->c[0].imag = c0i;
  c->c[1].real = c1r;
  c->c[1].imag = c1i;
  c->c[2].real = c2r;
  c->c[2].imag = c2i;

}

static inline void mult_adj_su3_mat_4vec( su3_matrix *mat, su3_vector *src,
			    su3_vector *dest0, su3_vector *dest1, 
			    su3_vector *dest2, su3_vector *dest3  ){
  register int n;
  register float c0r,c0i,c1r,c1i,c2r,c2i;
  register float br,bi,a0,a1,a2;
  register su3_matrix *a;
  register su3_vector *b,*c;
  su3_vector *cc[4] ;
  
  cc[0] = dest0 ; cc[1] = dest1 ;
  cc[2] = dest2 ; cc[3] = dest3 ;

  a = mat ; c=dest0 ; b = src;
  for(n=0;n<4;n++,a++,c=cc[n]){

  br=b->c[0].real;    bi=b->c[0].imag;
  a0=a->e[0][0].real;
  a1=a->e[0][1].real;
  a2=a->e[0][2].real;

  c0r = a0*br;
  c1r = a1*br;
  c2r = a2*br;
  c0i = a0*bi;
  c1i = a1*bi;
  c2i = a2*bi;

  a0=a->e[0][0].imag;
  a1=a->e[0][1].imag;
  a2=a->e[0][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[1].real;    bi=b->c[1].imag;
  a0=a->e[1][0].real;
  a1=a->e[1][1].real;
  a2=a->e[1][2].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=a->e[1][0].imag;
  a1=a->e[1][1].imag;
  a2=a->e[1][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  br=b->c[2].real;    bi=b->c[2].imag;
  a0=a->e[2][0].real;
  a1=a->e[2][1].real;
  a2=a->e[2][2].real;

  c0r += a0*br;
  c1r += a1*br;
  c2r += a2*br;
  c0i += a0*bi;
  c1i += a1*bi;
  c2i += a2*bi;

  a0=a->e[2][0].imag;
  a1=a->e[2][1].imag;
  a2=a->e[2][2].imag;

  c0r += a0*bi;
  c1r += a1*bi;
  c2r += a2*bi;
  c0i -= a0*br;
  c1i -= a1*br;
  c2i -= a2*br;

  c->c[0].real = c0r;
  c->c[0].imag = c0i;
  c->c[1].real = c1r;
  c->c[1].imag = c1i;
  c->c[2].real = c2r;
  c->c[2].imag = c2i;
  }
}
#endif
