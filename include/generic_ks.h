#ifndef _GENERIC_KS_H
#define _GENERIC_KS_H
/************************ generic_ks.h **********************************
*									*
*  Macros and declarations for generic_ks routines                      *
*  This header is for codes that call generic_ks routines               *
*  MIMD version 6 							*
*									*
*/

#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic_quark_types.h"
#include "../include/comdefs.h"

int congrad( int niter, float rsqmin, int parity, float *rsq );
void copy_latvec(field_offset src, field_offset dest, int parity);
void dslash( field_offset src, field_offset dest, int parity );
void dslash_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void clear_latvec(field_offset v,int parity);

void scalar_mult_latvec(field_offset src, float scalar,
			field_offset dest, int parity);
void scalar_mult_add_latvec(field_offset src1, field_offset src2,
			    float scalar, field_offset dest, int parity);
void scalar2_mult_add_su3_vector(su3_vector *a, float s1, su3_vector *b, 
				 float s2, su3_vector *c);

void scalar2_mult_add_latvec(field_offset src1,float scalar1,
			     field_offset src2,float scalar2,
			     field_offset dest,int parity);
void checkmul();
void phaseset();
void rephase( int flag );

void prefetch_vector( su3_vector * );
void prefetch_matrix( su3_matrix * );

int ks_congrad( field_offset src, field_offset dest, float mass,
     int niter, float rsqmin, int parity, float *rsq );

void cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]);
void cleanup_dslash_temps();
void dslash_fn( field_offset src, field_offset dest, int parity );
void ddslash_fn_du0( field_offset src, field_offset dest, int parity );
void dslash_fn_alltemp_special(su3_vector *src, su3_vector *dest,
			       int parity, msg_tag **tag, int start );
void dslash_fn_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );
void dslash_fn_on_temp( su3_vector *src, su3_vector *dest, int parity );
void ddslash_fn_du0_on_temp( su3_vector *src, su3_vector *dest, int parity );
void dslash_fn_on_temp_special(su3_vector *src, su3_vector *dest,
			       int parity, msg_tag **tag, int start );

void dslash_eo( field_offset src, field_offset dest, int parity );
void dslash_eo_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start );

int congrad_ks(            /* Return value is number of iterations taken */
     field_offset src,       /* type su3_vector* (preloaded source) */
     field_offset dest,      /* type su3_vector*  (answer and initial guess) */
     quark_invert_control *qic, /* inverter control */
     void *dmp               /* parameters defining the Dirac matrix */
     );

int ks_invert( /* Return value is number of iterations taken */
    field_offset src,   /* type su3_vector or multi_su3_vector 
			   (preloaded source) */
    field_offset dest,  /* type su3_vector or multi_su3_vector 
			   (answer and initial guess) */
    int (*invert_func)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    );

int ks_multicg(	/* Return value is number of iterations taken */
    field_offset src,	/* source vector (type su3_vector) */
    su3_vector **psim,	/* solution vectors */
    float *masses,	/* the masses */
    int num_masses,	/* number of masses */
    int niter,		/* maximal number of CG interations */
    float rsqmin,	/* desired residue squared */
    int parity,		/* parity to be worked on */
    float *final_rsq_ptr	/* final residue squared */
    );

/* f_meas.c */
void f_meas_imp( field_offset phi_off, field_offset xxx_off, float mass );

/* fpi_2.c */
int fpi_2( /* Return value is number of C.G. iterations taken */
  float *masses,   /* array of masses */
  int nmasses,      /* number of masses */
  float tol        /* tolerance for inverter check. */
  );

/* flavor_ops.c */
void sym_shift(int dir, field_offset src,field_offset dest) ;
void zeta_shift(int n, int *d, field_offset src, field_offset dest ) ;
void eta_shift(int n, int *d, field_offset src, field_offset dest ) ;


void mult_flavor_vector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_flavor_pseudoscalar(field_offset src, field_offset dest ) ;

void mult_spin_vector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest ) ;
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest ) ;
void mult_spin_pseudoscalar(field_offset src, field_offset dest ) ;

/* grsource.c */
void grsource(int parity);

/* grsource_imp.c */
void grsource_imp( field_offset dest, float mass, int parity);

/* mat_invert.c */
int mat_invert_cg( field_offset src, field_offset dest, field_offset temp,
		   float mass );
int mat_invert_uml(field_offset src, field_offset dest, field_offset temp,
		   float mass );
void check_invert( field_offset src, field_offset dest, float mass,
		   float tol);
/* nl_spectrum.c */
int nl_spectrum( float vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2);

/* quark_stuff.c */
void make_path_table();
void eo_fermion_force( float eps, int nflavors, field_offset x_off );
void eo_fermion_force_3f( float eps, int nflav1, field_offset x1_off,
	int nflav2, field_offset x2_off  );
void load_longlinks();
void load_fatlinks();

/* spectrum.c */
int spectrum();

/* spectrum2.c */
int spectrum2( float vmass, field_offset temp1, field_offset temp2 );

/* spectrum_hybrids.c */
int spectrum_hybrids( float mass, field_offset temp, float tol );

/* spectrum_mom.c */
int spectrum_mom( float qmass, float amass, field_offset temp, float tol);

/* spectrum_multimom.c */
int spectrum_multimom( float dyn_mass, float low_mass, float mass_inc, int nmasses, float tol);

/* spectrum_nd.c */
int spectrum_nd( float mass1, float mass2, float tol );

/* spectrum_nlpi2.c */
int spectrum_nlpi2( float qmass, float amass, field_offset temp, float tol);
void mult_rho0( int fdir,  field_offset src, field_offset dest ) ;
void mult_rhos( int fdir,  field_offset src, field_offset dest ) ;

#endif /* _GENERIC_KS_H */
