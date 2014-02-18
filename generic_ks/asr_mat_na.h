inline void asr_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
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

void asr_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
  int i,j;
  register float t,ar,ai,br,bi,cr,ci;
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
}

void asr_su3_an( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
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
