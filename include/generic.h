#ifndef _GENERIC_H
#define _GENERIC_H
/************************ generic.h *************************************
*									*
*  Macros and declarations for miscellaneous generic routines           *
*  This header is for codes that call generic routines                  *
*  MIMD version 6 							*
*									*
*/

/* Other generic directory declarations are elsewhere:

   For com_*.c, see comdefs.h
   For io_lat4.c io_ansi.c, io_nonansi.c, io_piofs.c, io_romio.c see io_lat.h
   For io_wb3.c, see io_wb.h
*/

#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/macros.h"
#include "../include/random.h"

/* ape_smear.c */
void ape_smear_dir(
  field_offset src,       /* field offset for su3_matrix[4] type 
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for su3_matrix type 
			     pointing to a specific direction 
			     output smeared links */
  float staple_weight,    /* single staple weight */
  float link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  float tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  );

void ape_smear(
  field_offset src,       /* field offset for su3_matrix type 
			     input unsmeared links */
  field_offset dest,      /* field offset for su3_matrix type 
			     output smeared links */
  float staple_weight,    /* single staple weight */
  float link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
 			          only spacelike staples 
			     = 0 (false) smear all links with
			     all staples */
  int nhits,              /* reproject onto SU(3): number of 
			     SU(2) hits. 0 for no reprojection */
  float tol               /* tolerance for SU(3) projection.
			     If nonzero, treat nhits as a maximum
			     number of hits.  If zero, treat nhits
			     as a prescribed number of hits. */ 
  );

/* ax_gauge.c */
void ax_gauge();

/* bsd_sum.c */
int32type bsd_sum (char *data,int32type total_bytes);

/* check_unitarity.c */
float check_unitarity( void );

/* d_plaq?.c */
void d_plaquette(double *ss_plaq,double *st_plaq);

/* field_strength.c */
void make_field_strength(
  field_offset link_src,       /* field offset for su3_matrix[4] type 
				  for the source link matrices */
  field_offset field_dest      /* field offset for su3_matrix[6] type
				  for the resulting field strength */
  );

/* gaugefix.c and gaugefix2.c */
void gaugefix(int gauge_dir,float relax_boost,int max_gauge_iter,
	      float gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] );

/* gauge_stuff.c */
double imp_gauge_action();
void imp_gauge_force( float eps, field_offset mom_off );
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);

/* glueball_op.c */
void make_glueball_ops();
void measure_glueball_ops();

/* hvy_pot.c */
void hvy_pot( field_offset links );

/* layout_*.c */
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
int num_sites(int node);

/* make_lattice.c */
void make_lattice();

/* make_global_fields.c */
void make_global_fields();

/* path_product.c */
void path_product( const int *dir, const int length);
void path_prod_subl(const int *dir, const int length, const int subl);

/* plaquette4.c */
void plaquette(float *ss_plaq,float *st_plaq);

/* ploop?.c */
complex ploop( void );

/* ploop_staple.c */
complex ploop_staple(float alpha_fuzz);

/* project_su3_hit.c */
void project_su3(
   su3_matrix *w,         /* input initial guess. output resulting
                             SU(3) matrix */
   su3_matrix *q,         /* starting 3 x 3 complex matrix */
   int Nhit,              /* number of SU(2) hits. 0 for no projection */
   float tol              /* tolerance for SU(3) projection.
			     If nonzero, treat Nhit as a maximum
			     number of hits.  If zero, treat Nhit
			     as a prescribed number of hits. */ 
   );

/* rand_gauge.c */
void rand_gauge(field_offset G);

/* ranmom.c */
void ranmom( void );

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
float myrand(double_prn *prn_pt);

/* restrict_fourier.c */
void setup_restrict_fourier( int *key, int *slice);
void restrict_fourier( 
     field_offset src,	 /* src is field to be transformed */
     field_offset space, /* space is working space, same size as src */
     field_offset space2,/* space2 is working space, same size as src */
                         /* space2 is needed only for non power of 2 */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign);	 /* 1 for x -> k, -1 for k -> x */

/* reunitarize2.c */
void reunitarize( void );
int reunit_su3(su3_matrix *c);

#endif	/* _GENERIC_H */
