static INLINE void mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c ){
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

