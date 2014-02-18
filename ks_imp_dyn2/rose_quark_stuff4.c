/****** quark_stuff.c  -- ******************/
/* MIMD version 6 */
/* quark action stuff for improved action
* D.T. 1/28/98, starting from gauge_stuff.c
* K.O. 3/99 Added optimized fattening for Asq actions
* D.T. 4/99 Combine force calculations for both mass quarks
* K.O. 4/99 Optimized force for Asq action
* S.G. 7/01, modified to use t_longlink and t_fatlink
* C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
* T.B. 11/01, Added d(M)/d(u0) dependencies for equation of state calc's with
*             Asqtad action - ATTN: site structure needs 'dfatlink_du0' in
*             addition to 'fatlink': #define DM_DU0
*
* J.O. 3/04 Rearranged loops for optimization
* J.O. C.D. 3/04 Copied forward links for optimization and 
*                kept mtags open for restart_gather
*                Worked with pointers where possible to avoid copying.
* This code combines quark_stuff.c and quark_stuff_tmp.c
* with DSLASH_TMP_LINKS defined, puts links in field-major variables
* t_longlink and t_fatlink. Otherwise, puts them in the site structure.
*
* In this directory, assume all paths connect even to odd sites, etc.
* Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
* as separate parity transforms of the fundamental paths.  They will
* generally need a negative sign in Dslash.  See bottom for a long
* comment on sign conventions.
*/
/*
* 9/10/01, flopcount for ASQ_OPTIMIZED
* Fermion force: 420528 for eo_fermion_force_3f() (need to check this)
* Load fatlinks: 50616
*/
/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fermion force: 253935 for eo_fermion_force()
 * Fermion force: 433968 for eo_fermion_force_3f()
 * Fatlinks:       61632 for load_fatlinks
 */
/*#define FFTIME*/
/*#define LLTIME*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "asr_mat_na.h"
#define NULL_FP -1 /* NULL field_offset to be used in the optimized version *
                    * of the load_fatlinks subroutine */
/* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom. */
#include <quark_action.h>
/* Include file specifies the basic paths */
#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void printpath(int *path,int length);
void path_transport(field_offset src,field_offset dest,int parity,int *dir,int length);
void path_transport_hwv(field_offset src,field_offset dest,int parity,int *dir,int length);
#ifdef  ASQ_OPTIMIZED_FATTENING
#ifdef DM_DU0
#ifndef TADPOLE_IMPROVE
#endif
#else
void compute_gen_staple(field_offset staple,int mu,int nu,field_offset link,float coef);
#endif
#endif
#ifdef  ASQ_OPTIMIZED_FORCE
#ifndef FN
#endif
void u_shift_fermion(su3_vector *src,su3_vector *dest,int dir);
void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,float coef);
void side_link_force(int mu,int nu,float coeff,su3_vector *Path,su3_vector *Path_nu,su3_vector *Path_mu,su3_vector *Path_numu);
void u_shift_hw_fermion_np(half_wilson_vector *src,half_wilson_vector **dest_pt,int dir,msg_tag **mtag,half_wilson_vector *tmpvec);
void u_shift_hw_fermion_pp(half_wilson_vector **src_pt,half_wilson_vector **dest_pt,int dir,msg_tag **mtag,half_wilson_vector *tmpvec);
void add_3f_force_to_mom_nn(half_wilson_vector *back,half_wilson_vector *forw,int dir,float coeff[2UL]);
void add_3f_force_to_mom_np(half_wilson_vector *back,half_wilson_vector **forw_pt,int dir,float (coeff)[2UL]);
void add_3f_force_to_mom_pn(half_wilson_vector **back_pt,half_wilson_vector *forw,int dir,float (coeff)[2UL]);
void add_3f_force_to_mom_pp(half_wilson_vector **back_pt,half_wilson_vector **forw_pt,int dir,float (coeff)[2UL]);
void side_link_3f_force_pppp(int mu,int nu,float coeff[2UL],half_wilson_vector **Path,half_wilson_vector **Path_nu,half_wilson_vector **Path_mu,half_wilson_vector **Path_numu);
void side_link_3f_force_npnp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector **Path_nu,half_wilson_vector *Path_mu,half_wilson_vector **Path_numu);
void side_link_3f_force_nnpp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector *Path_nu,half_wilson_vector **Path_mu,half_wilson_vector **Path_numu);
void side_link_3f_force_nnnp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector *Path_nu,half_wilson_vector *Path_mu,half_wilson_vector **Path_numu);
#endif
/* number of rotations/reflections for each 
					kind */
int path_num[6UL];
/* actual path coefficient     *
                                               * it is equal to path_coeff   *
                                               * if not tadpole improvement  *
                                               * is specified                *
                                               * or path_coeff*u_0^(L-1) when*
                                               * tadpole improvement is      *
                                               * specified                   */
static float act_path_coeff[6UL];
#ifdef DM_DU0
/* coefficient for
						  d(Dslash)/d(u0) */
#endif
/* Array of structures, for each rotation and reflection of each kind of
	path.  */
struct __unnamed_class___F0_L144_C1_unknown_scope_and_name_variable_declaration__variable_type__Ab_i_index_7_Ae__variable_name_unknown_scope_and_name__scope__dir__DELIMITER__unknown_scope_and_name_variable_declaration__variable_type_i_variable_name_unknown_scope_and_name__scope__length__DELIMITER__unknown_scope_and_name_variable_declaration__variable_type_f_variable_name_unknown_scope_and_name__scope__coeff__DELIMITER__unknown_scope_and_name_variable_declaration__variable_type_f_variable_name_unknown_scope_and_name__scope__forwback {
/* directions in path */
int dir[7UL];
/* length of path */
int length;
/* coefficient, including minus sign if backwards */
float coeff;
#ifdef DM_DU0
/* coefficient for d(Dslash)/d(u0) */
#endif
/* +1 if in forward Dslash, -1 if in backward */
float forwback;}q_paths[688UL];
/* number of paths in dslash */
int num_q_paths;
/* number of paths before rotation/reflection */
int num_basic_paths;
int is_path_equal(int *path1,int *path2,int length);
int add_basic_path(int *vec,int length,float coeff);
#ifdef QSINLINE
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) {\
\
  float _a0r,_a0i,_a1r,_a1i,_a2r,_a2i;\
  float _b0r,_b0i,_b1r,_b1i,_b2r,_b2i;\
  \
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
}
#define su3_projector_for_inline( a, b, dest ){\
  int _i,_j;\
  float _tmp,_tmp2;\
    for(_i=0;_i<3;_i++)for(_j=0;_j<3;_j++){\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].real;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].imag;\
	(dest)->e[_i][_j].real = _tmp + _tmp2;\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].imag;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].real;\
	(dest)->e[_i][_j].imag = _tmp - _tmp2;\
    }\
}
#else /* External versions */
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) mult_su3_mat_hwvec( mat, src, dest ) 
#define su3_projector_for_inline( a, b, dest ) su3_projector( a, b, dest )
#endif
/* Make table of paths in action */

void make_path_table()
{
  int i;
  int j;
#ifdef TADPOLE_IMPROVE
  int k;
#endif
/* table of directions, 1 for each kind of path */
/**int path_ind[MAX_BASIC_PATHS][MAX_LENGTH];**/
/* table of coefficients in action, for each path */
  if (this_node == 0) 
    printf("%s\n",quark_action_description);
  num_q_paths = 0;
  num_basic_paths = 0;
/* add rots. and reflects to table, print out the path coefficients */
  if (this_node == 0) 
    printf("path coefficients: npath  path_coeff  multiplicity\n");
  for (j = 0; j < quark_action_npaths; j++) {
    float this_coeff;
    this_coeff = path_coeff[j];
#ifdef TADPOLE_IMPROVE
    for (k = 1; k < path_length_in[j]; k++) 
      this_coeff /= u0;
#endif
    act_path_coeff[j] = this_coeff;
#ifdef DM_DU0
#endif
    i = add_basic_path(path_ind[j],path_length_in[j],this_coeff);
    if (this_node == 0) 
      printf("                    %d      %e     %d\n",j,this_coeff,i);
  }
}
/* add rotations and reflections of a path to the table.  Return
   multiplicity of paths added */

int add_basic_path(int *basic_vec,int length,float coeff)
{
  int perm[8UL];
  int pp[8UL];
  int ir[4UL];
  int j;
  int path_num;
  int vec[7UL];
  int flag;
  path_num = 0;
/* now fill the long table with all rotations and reflections
	of the fundamental path.  The path presented to us is for
        the positive x component of dslash, so if the x coordinate
        is reflected it will appear with a negative sign. */
/* permutations */
  for (perm[0] = 0; perm[0] < 4; perm[0]++) 
    for (perm[1] = 0; perm[1] < 4; perm[1]++) 
      for (perm[2] = 0; perm[2] < 4; perm[2]++) 
        for (perm[3] = 0; perm[3] < 4; perm[3]++) {
          if ((((((perm[0] != perm[1]) && (perm[0] != perm[2])) && (perm[0] != perm[3])) && (perm[1] != perm[2])) && (perm[1] != perm[3])) && (perm[2] != perm[3])) {
/* reflections*/
            for (ir[0] = 0; ir[0] < 2; ir[0]++) 
              for (ir[1] = 0; ir[1] < 2; ir[1]++) 
                for (ir[2] = 0; ir[2] < 2; ir[2]++) 
                  for (ir[3] = 0; ir[3] < 2; ir[3]++) {
                    for (j = 0; j < 4; j++) {
                      pp[j] = perm[j];
                      if (ir[j] == 1) 
                        pp[j] = (7 - pp[j]);
                      pp[7 - j] = (7 - pp[j]);
                    }
/* create new vector*/
                    for (j = 0; j < length; j++) 
                      vec[j] = pp[basic_vec[j]];
                    for (j = length; j < 7; j++) 
                      vec[j] = -1;
                    flag = 0;
{
/* check if it's a new set: */
                      for (j = 0; j < num_q_paths; j++) {
                        flag = is_path_equal(vec,q_paths[j].dir,7);
                        if (flag == 1) 
                          break; 
                      }
                    }
                    if (flag == 0) {
                      if (num_q_paths >= 688) {
                        if (this_node == 0) 
                          printf("OOPS: MAX_NUM too small\n");
                        exit(0);
                      }
                      q_paths[num_q_paths].length = length;
                      for (j = 0; j < 7; j++) 
                        q_paths[num_q_paths].dir[j] = vec[j];
/* remember to copy NODIR's, or comparison will get confused */
                      if (ir[0] == 0) {
                        q_paths[num_q_paths].coeff = coeff;
#ifdef DM_DU0
#endif
                        q_paths[num_q_paths].forwback = (+1);
                      }
                      else {
                        q_paths[num_q_paths].coeff = -coeff;
#ifdef DM_DU0
#endif
                        q_paths[num_q_paths].forwback = (-1);
                      }
                      num_q_paths++;
                      path_num++;
/**node0_printf("ADD PATH %d:  rx=%d ",num_q_paths-1,ir[0]);
		 printpath( vec, length );**/
                    }
/* end reflection*/
                  }
/* end permutation if block */
          }
/* end permutation */
        }
  num_basic_paths++;
  return path_num;
/* add_basic_path */
}
/* parallel transport a vector to the current site along a path.
   For example, if the path is "XUP", bring in the vector from
   the +X direction to the current site. If KS phases are in lattice,
   this transport will automatically include them. 
   OK for src and dest to be the same.  OK for length=0.  */
/* NOT OPTIMIZED at the moment - do lots of extra copying.  Use temp.
   vectors rather than stuff in lattice.h */

void path_transport(field_offset src,field_offset dest,int parity,int *dir,int length)
{
  register int i;
  register site *s;
  msg_tag *mtag0;
  int j;
/*source, dest and workspace*/
  su3_vector *tmp_src;
  su3_vector *tmp_dest;
  su3_vector *tmp_work;
/* scratch */
  su3_vector *tmp_pt;
/* parity for this step */
  int tmp_parity;
  int tmp_otherparity;
  if (length > 0) {
    tmp_src = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    tmp_dest = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    tmp_work = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    for (j = (length - 1); j >= 0; j--) {
/* figure out parities for this step */
      if ((j % 2) == 0) {
        tmp_parity = parity;
        switch(tmp_parity){
          case 0x02:
{
            tmp_otherparity = 1;
            break; 
          }
          case 0x01:
{
            tmp_otherparity = 2;
            break; 
          }
          case 0x03:
{
            tmp_otherparity = 3;
            break; 
          }
        }
      }
      else 
/* odd # step */
{
        tmp_otherparity = parity;
        switch(tmp_otherparity){
          case 0x02:
{
            tmp_parity = 1;
            break; 
          }
          case 0x01:
{
            tmp_parity = 2;
            break; 
          }
          case 0x03:
{
            tmp_parity = 3;
            break; 
          }
        }
      }
      if (j == (length - 1)) {
        for (((i = ((tmp_otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_otherparity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          tmp_src[i] =  *((su3_vector *)(((char *)s) + src));
        }
      }
      if (dir[j] <= 3) {
        mtag0 = start_gather_from_temp(tmp_src,(sizeof(su3_vector )),dir[j],tmp_parity,gen_pt[0]);
        wait_gather(mtag0);
        for (((i = ((tmp_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          mult_su3_mat_vec(((s -> link) + dir[j]),((su3_vector *)gen_pt[0][i]),(tmp_dest + i));
        }
        cleanup_gather(mtag0);
      }
      else 
/* GOES_BACKWARDS(dir[j]) */
{
        for (((i = ((tmp_otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_otherparity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          mult_adj_su3_mat_vec(((s -> link) + (7 - dir[j])),(tmp_src + i),(tmp_work + i));
        }
        mtag0 = start_gather_from_temp(tmp_work,(sizeof(su3_vector )),dir[j],tmp_parity,gen_pt[0]);
        wait_gather(mtag0);
        for (((i = ((tmp_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          tmp_dest[i] =  *((su3_vector *)gen_pt[0][i]);
        }
        cleanup_gather(mtag0);
      }
/* src for next step is dest for this one. */
      tmp_pt = tmp_src;
      tmp_src = tmp_dest;
      tmp_dest = tmp_pt;
/* j=link in path */
    }
/* done, copy result into real dest. (tmp_src now points to result) */
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
       *((su3_vector *)(((char *)s) + dest)) = tmp_src[i];
    }
    free(tmp_src);
    free(tmp_dest);
    free(tmp_work);
/* end if(length>0) */
  }
  else 
/* for length=0 */
if (src != dest) {
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
       *((su3_vector *)(((char *)s) + dest)) =  *((su3_vector *)(((char *)s) + src));
    }
  }
/* path_transport */
}
/* Path transport a half_wilson_vector */

void path_transport_hwv(field_offset src,field_offset dest,int parity,int *dir,int length)
{
  register int i;
  register site *s;
  msg_tag *mtag0;
  int j;
/*source, dest and workspace*/
  half_wilson_vector *tmp_src;
  half_wilson_vector *tmp_dest;
  half_wilson_vector *tmp_work;
/* scratch */
  half_wilson_vector *tmp_pt;
/* parity for this step */
  int tmp_parity;
  int tmp_otherparity;
  if (length > 0) {
    tmp_src = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    tmp_dest = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    tmp_work = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    for (j = (length - 1); j >= 0; j--) {
/* figure out parities for this step */
      if ((j % 2) == 0) {
        tmp_parity = parity;
        switch(tmp_parity){
          case 0x02:
{
            tmp_otherparity = 1;
            break; 
          }
          case 0x01:
{
            tmp_otherparity = 2;
            break; 
          }
          case 0x03:
{
            tmp_otherparity = 3;
            break; 
          }
        }
      }
      else 
/* odd # step */
{
        tmp_otherparity = parity;
        switch(tmp_otherparity){
          case 0x02:
{
            tmp_parity = 1;
            break; 
          }
          case 0x01:
{
            tmp_parity = 2;
            break; 
          }
          case 0x03:
{
            tmp_parity = 3;
            break; 
          }
        }
      }
      if (j == (length - 1)) {
        for (((i = ((tmp_otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_otherparity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          tmp_src[i] =  *((half_wilson_vector *)(((char *)s) + src));
        }
      }
      if (dir[j] <= 3) {
        mtag0 = start_gather_from_temp(tmp_src,(sizeof(half_wilson_vector )),dir[j],tmp_parity,gen_pt[0]);
        wait_gather(mtag0);
        for (((i = ((tmp_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          mult_su3_mat_hwvec(((s -> link) + dir[j]),((half_wilson_vector *)gen_pt[0][i]),(tmp_dest + i));
        }
        cleanup_gather(mtag0);
      }
      else 
/* GOES_BACKWARDS(dir[j]) */
{
        for (((i = ((tmp_otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_otherparity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          mult_adj_su3_mat_hwvec(((s -> link) + (7 - dir[j])),(tmp_src + i),(tmp_work + i));
        }
        mtag0 = start_gather_from_temp(tmp_work,(sizeof(half_wilson_vector )),dir[j],tmp_parity,gen_pt[0]);
        wait_gather(mtag0);
        for (((i = ((tmp_parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((tmp_parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
          tmp_dest[i] =  *((half_wilson_vector *)gen_pt[0][i]);
        }
        cleanup_gather(mtag0);
      }
/* src for next step is dest for this one. */
      tmp_pt = tmp_src;
      tmp_src = tmp_dest;
      tmp_dest = tmp_pt;
/* j=link in path */
    }
/* done, copy result into real dest. (tmp_src now points to result) */
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
       *((half_wilson_vector *)(((char *)s) + dest)) = tmp_src[i];
    }
    free(tmp_src);
    free(tmp_dest);
    free(tmp_work);
/* end if(length>0) */
  }
  else 
/* for length=0 */
if (src != dest) {
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < (((parity == 2)?even_sites_on_node : sites_on_node)); (i++ , s++)) {
       *((half_wilson_vector *)(((char *)s) + dest)) =  *((half_wilson_vector *)(((char *)s) + src));
    }
  }
/* path_transport_hwv */
}
#ifdef EO
/* Stupid dslash routine that follows all the paths.
  Should optimize by precomputing sums of paths to all the displacements
  Use dslash_fn for actions that involve only +X and +X+X+X couplings.
*/
/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  */
/* coefficient of path */
/* Parallel transport by all the paths in the action.  
       Multiply by coefficient in table
    */
/* loop over paths */
/* ipath */
/* dslash_eo */
#endif /* ifdef EO */
#ifdef FN
/* Sum over paths connecting to nearest neighbor point (fat link) and to third
   nearest neighbor (longlinks) */
/* Doug Toussaint 2/4/98 */
/* modified to use t_longlinks, S. Gottlieb 7/13/01 */
/* long link calculating routine */
/* path_product() follows the path starting at step 0, and
   leaves the answer at the end of the path.  We want the answer
   at the site where the path begins.  So we look for paths with
   the opposite displacement from the displacement of the point
   that we want to transport to this site, and take the adjoint
   of the matrix at the end. clear? */
/* KS phases and APBC must be in the links. See long comment at bottom*/

void load_longlinks()
{
  register int i;
  register site *s;
  int ipath;
  int dir;
  int disp[4UL];
  int nflop = 1804;
  register su3_matrix *long1;
#ifdef LLTIME
  double dtime;
  dtime = -dclock();
#endif
  if (phases_in != 1) {
    if (this_node == 0) 
      printf("BOTCH: load_longlinks needs phases in\n");
    terminate(0);
  }
/* loop over longlink directions */
  for (dir = 0; dir <= 3; dir++) {
/* set longlink to zero */
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
#ifdef DSLASH_TMP_LINKS
      long1 = (t_longlink + ((4 * i) + dir));
#else
#endif
      clear_su3mat(long1);
    }
/* loop over paths, checking for ones with total displacement 3*dir */
/* loop over paths */
    for (ipath = 0; ipath < num_q_paths; ipath++) {
/* compute total displacement of path */
      for (i = 0; i <= 3; i++) 
        disp[i] = 0;
      for (i = 0; i < q_paths[ipath].length; i++) {
        if (q_paths[ipath].dir[i] <= 3) 
          disp[q_paths[ipath].dir[i]]++;
        else 
          disp[7 - q_paths[ipath].dir[i]]--;
      }
{
        for (((disp[dir] += 3) , (i = 0)); i <= 3; i++) 
          if (disp[i] != 0) 
            break; 
      }
/* skip if path doesn't go to right place */
      if (i <= 3) 
        continue; 
/**printf("ipath = %d, found a path:  ",ipath);
for(j=0;j<q_paths[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/
      path_product(q_paths[ipath].dir,q_paths[ipath].length);
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        su3_adjoint(&s -> tempmat1,&s -> staple);
#ifdef DSLASH_TMP_LINKS
        long1 = (t_longlink + ((4 * i) + dir));
#else
#endif
        scalar_mult_add_su3_matrix(long1,&s -> staple,-q_paths[ipath].coeff,long1);
/* minus sign in coeff. because we used backward path*/
      }
/* ipath */
    }
/* loop over directions */
  }
  valid_longlinks = 1;
#ifdef LLTIME
  dtime += dclock();
  if (this_node == 0) 
    printf("LLTIME(long): time =  %e (Naik) mflops = %e\n",dtime,((((float )nflop) * volume) / ((1e6 * dtime) * (numnodes()))));
#endif
/* load_longlinks() */
}
/* KS phases and APBC must be in the links. See long comment at bottom*/

void load_fatlinks()
{
  register int i;
  register site *s;
  int dir;
  register su3_matrix *fat1;
#ifdef DM_DU0
#endif
#ifdef ASQ_OPTIMIZED_FATTENING
  int nu;
  int rho;
  int sig;
/* needed to fix the problem with the Lepage
			       term */
  float one_link;
  float one_link2;
#else
#endif
#ifdef LLTIME
  int nflop = 61632;
  double dtime;
  dtime = -dclock();
#endif
  if (phases_in != 1) {
    if (this_node == 0) 
      printf("BOTCH: load_fatlinks needs phases in\n");
    terminate(0);
  }
#ifndef  ASQ_OPTIMIZED_FATTENING   /* general case code */
/* loop over fatlink directions */
/* set fatlink to zero */
#ifdef DSLASH_TMP_LINKS
#else
#endif
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
#endif
/* loop over paths, checking for ones with total displacement 1*dir */
/* loop over paths */
/* compute total displacement of path */
/* skip if path doesn't go to right place */
/**printf("dir = %d, found a path:  ",dir);
for(j=0;j<q_paths.[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* minus sign in coeff. because we used backward path*/
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
/* minus sign in coeff. because we used backward path*/
#endif
/* ipath */
/* loop over directions */
#else	/* ASQ_OPTIMIZED_FATTENING, for Asq and Asqtad actions */
/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
/* to fix up the Lepage term, included by a trick below */
  one_link = (act_path_coeff[0] - (6.0 * act_path_coeff[5]));
#ifdef DM_DU0
#endif
  for (dir = 0; dir <= 3; dir++) {
/* Intialize fat links with c_1*U_\mu(x) */
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
#ifdef DSLASH_TMP_LINKS
      fat1 = (t_fatlink + ((4 * i) + dir));
#else
#endif
      scalar_mult_su3_matrix(((s -> link) + dir),one_link,fat1);
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
#endif
    }
    for (nu = 0; nu <= 3; nu++) 
      if (nu != dir) {
#ifdef DM_DU0
/* The Lepage term */
/* Note this also involves modifying c_1 (above) */
/* sig */
/* rho */
#else
        compute_gen_staple(((field_offset )(((char *)(&lattice[0].staple)) - ((char *)(lattice + 0)))),dir,nu,((field_offset )(((char *)(lattice[0].link + dir)) - ((char *)(lattice + 0)))),act_path_coeff[2]);
/* The Lepage term */
/* Note this also involves modifying c_1 (above) */
        compute_gen_staple((-1),dir,nu,((field_offset )(((char *)(&lattice[0].staple)) - ((char *)(lattice + 0)))),act_path_coeff[5]);
        for (rho = 0; rho <= 3; rho++) 
          if ((rho != dir) && (rho != nu)) {
            compute_gen_staple(((field_offset )(((char *)(&lattice[0].tempmat1)) - ((char *)(lattice + 0)))),dir,rho,((field_offset )(((char *)(&lattice[0].staple)) - ((char *)(lattice + 0)))),act_path_coeff[3]);
            for (sig = 0; sig <= 3; sig++) 
              if (((sig != dir) && (sig != nu)) && (sig != rho)) {
                compute_gen_staple((-1),dir,sig,((field_offset )(((char *)(&lattice[0].tempmat1)) - ((char *)(lattice + 0)))),act_path_coeff[4]);
/* sig */
              }
/* rho */
          }
#endif
/* nu */
      }
/* dir */
  }
#endif
  valid_fatlinks = 1;
#ifdef LLTIME
  dtime += dclock();
  if (this_node == 0) 
    printf("LLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,((((float )nflop) * volume) / ((1e6 * dtime) * (numnodes()))));
#endif
/* load_fatlinks() */
}
#endif /* ifdef FN */
/* compare two paths, return 1 if equal, else zero */

int is_path_equal(int *path1,int *path2,int length)
{
  register int i;
  for (i = 0; i < length; i++) 
    if (path1[i] != path2[i]) 
      return 0;
  return 1;
}
/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   x_off, and dslash(x_off,x_off,ODD) has been run. (fills in x_off_odd) */
/* SEE LONG COMMENTS AT END */
#ifndef ASQ_OPTIMIZED_FORCE
/* note CG_solution and Dslash * solution are combined in "x_off" */
/* New version 1/21/99.  Use forward part of Dslash to get force */
/* see long comment at end */
/* For each link we need x_off transported from both ends of path. */
#ifdef FFTIME
#endif
#ifdef FFTIME
#endif
/* Use half_wilson_vectors to store x_off transported from ends of
     path.  0 component from forward end, 1 component from back end */
/* loop over paths, and loop over links in path */
/* skip backwards dslash */
/* path transport x_off and Dslash*x_off from far end.  Sometimes
	we need them at the start point of the path, and sometimes
	one link into the path - an optimization for later */
/* use tempvec[1] for transport from starting end */
/* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
/* path transport x_off and Dslash*x_off from previous point */
/* Use "half_wilson_vector" to handle pair of vectors -
	0 component is x_off from forward end, 1 component from back end */
/* sometimes we don't need them */
/* GOES_BACKWARDS(lastdir) */
/* add in contribution to the force */
/* Put antihermitian traceless part into momentum */
/* end loop over links in path */
/* end loop over paths */
#ifdef FFTIME
/**printf("TLENGTH: %d\n",tlength);**/
#endif
/* eo_fermion_force(version 6) */
/* note CG_solution and Dslash * solution are combined in "x_off" */
/* New version 1/21/99.  Use forward part of Dslash to get force */
/* 4/14/99 combine force from two different mass quarks, (eg 2+1flavors) */
/* see long comment at end */
/* For each link we need x_off transported from both ends of path. */
#ifdef FFTIME
#endif
#ifdef FFTIME
#endif
/* Use wilson_vectors to store x_off transported from ends of
     path.  0 and 1 components from forward end, 2 and 3 components
     from back end */
/* loop over paths, and loop over links in path */
/* skip backwards dslash */
/* path transport x_off and Dslash*x_off from far end.  Sometimes
	we need them at the start point of the path, and sometimes
	one link into the path - an optimization for later */
/**
    path_transport( x1_off, F_OFFSET(tempvec[0]),
      EVENANDODD, q_paths[ipath].dir, length );
    path_transport( x2_off, F_OFFSET(tempvec[1]),
      EVENANDODD, q_paths[ipath].dir, length );
    **/
/** WARNING!! Assumes xxx1 and xxx2 contiguous **/
/* use tempvec[2] for transport from starting end */
/* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
/* path transport x_off and Dslash*x_off from previous point */
/* Use "wilson_vector" to handle pair of vectors -
	0 component is x_off1 from forward end, 1 component is x_off2
	from forward end,  2 and 3  components are x_off1 and x_off2 
	from back end */
/* sometimes we don't need them */
/* GOES_BACKWARDS(lastdir) */
/* add in contribution to the force */
/* Put antihermitian traceless part into momentum */
/* end loop over links in path */
/* end loop over paths */
#ifdef FFTIME
/**printf("TLENGTH: %d\n",tlength);**/
#endif
/* eo_fermion_force_3f(version 6) */
#else /* ASQ_OPTIMIZED_FORCE, for Asq and Asqtad actions */
/* Optimized force code for the Asq and Asqtad actions                 *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
#define Pmu          tempvec[0] 
#define Pnumu        tempvec[1]
#define Prhonumu     tempvec[2]
#define P7           tempvec[3]
#define P7rho        tempvec[4]              
#define P7rhonu      tempvec[5]
#define P5           tempvec[6]
#define P3           tempvec[7]
#define P5nu         tempvec[3]
#define P3mu         tempvec[3]
#define Popmu        tempvec[4]
#define Pmumumu      tempvec[4]

void eo_fermion_force(float eps,int nflavors,field_offset x_off)
{
/* note CG_solution and Dslash * solution are combined in "x_off" */
/* New version 1/21/99.  Use forward part of Dslash to get force */
/* see long comment at end */
/* For each link we need x_off transported from both ends of path. */
  register int i;
  register site *s;
  int mu;
  int nu;
  int rho;
  int sig;
  int DirectLinks[8UL];
  float ferm_epsilon;
  float coeff;
  float OneLink;
  float Lepage;
  float Naik;
  float FiveSt;
  float ThreeSt;
  float SevenSt;
  su3_vector *tempvec[8UL];
  su3_vector *temp_x;
#ifdef FFTIME
  int nflop = 253935;
  double dtime;
  dtime = -dclock();
#endif
  ferm_epsilon = ((2.0 * (nflavors / 4.0)) * eps);
/* Path coefficients times fermion epsilon */
  OneLink = (act_path_coeff[0] * ferm_epsilon);
  Naik = (act_path_coeff[1] * ferm_epsilon);
  ThreeSt = (act_path_coeff[2] * ferm_epsilon);
  FiveSt = (act_path_coeff[3] * ferm_epsilon);
  SevenSt = (act_path_coeff[4] * ferm_epsilon);
  Lepage = (act_path_coeff[5] * ferm_epsilon);
/* *************************************** */
/* Initialize the DirectLink flags */
  for (mu = 0; mu < 8; mu++) 
    DirectLinks[mu] = 0;
/* Allocate temporary vectors */
  for (mu = 0; mu < 8; mu++) 
    tempvec[mu] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
/*copy x_off to a temporary vector */
  temp_x = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) 
    temp_x[i] =  *((su3_vector *)(((char *)s) + x_off));
  for (sig = 0; sig < 8; sig++) {
    for (mu = 0; mu < 8; mu++) 
      if ((mu != sig) && (mu != (7 - sig))) {
        u_shift_fermion(temp_x,tempvec[0],(7 - mu));
        u_shift_fermion(tempvec[0],tempvec[7],sig);
        if (sig <= 3) {
/* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
          add_force_to_mom(tempvec[7],tempvec[0],sig,-ThreeSt);
        }
        for (nu = 0; nu < 8; nu++) 
          if ((((nu != mu) && (nu != (7 - mu))) && (nu != sig)) && (nu != (7 - sig))) {
            u_shift_fermion(tempvec[0],tempvec[1],(7 - nu));
            u_shift_fermion(tempvec[1],tempvec[6],sig);
            if (sig <= 3) {
/* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
              add_force_to_mom(tempvec[6],tempvec[1],sig,FiveSt);
            }
            for (rho = 0; rho < 8; rho++) 
              if ((((((rho != mu) && (rho != (7 - mu))) && (rho != nu)) && (rho != (7 - nu))) && (rho != sig)) && (rho != (7 - sig))) {
                u_shift_fermion(tempvec[1],tempvec[2],(7 - rho));
/* Length 7 paths */
                u_shift_fermion(tempvec[2],tempvec[3],sig);
                if (sig <= 3) {
/* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
                  add_force_to_mom(tempvec[3],tempvec[2],sig,-SevenSt);
                }
/*Add the force F_rho the 2(4) link in the path: +     */
                u_shift_fermion(tempvec[3],tempvec[4],rho);
                side_link_force(rho,sig,SevenSt,tempvec[1],tempvec[3],tempvec[2],tempvec[4]);
/* Add the P7rho vector to P5 */
                coeff = (SevenSt / FiveSt);
                for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) 
                  scalar_mult_add_su3_vector((tempvec[6] + i),(tempvec[4] + i),coeff,(tempvec[6] + i));
/* rho */
              }
/* Length 5 paths */
/*Add the force F_nu the 1(3) link in the path: -     */
            u_shift_fermion(tempvec[6],tempvec[3],nu);
            side_link_force(nu,sig,-FiveSt,tempvec[0],tempvec[6],tempvec[1],tempvec[3]);
/* Add the P5nu vector to P3 */
            coeff = (FiveSt / ThreeSt);
            for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) 
              scalar_mult_add_su3_vector((tempvec[7] + i),(tempvec[3] + i),coeff,(tempvec[7] + i));
/* nu */
          }
/* Now the Lepage term... It is the same with 5-link paths with
             nu=mu and FiveSt=Lepage. So Pnumu is really Pmumu */
        u_shift_fermion(tempvec[0],tempvec[1],(7 - mu));
        u_shift_fermion(tempvec[1],tempvec[6],sig);
        if (sig <= 3) {
/* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
          add_force_to_mom(tempvec[6],tempvec[1],sig,Lepage);
        }
/*Add the force F_nu the 1(3) link in the path: -     */
        u_shift_fermion(tempvec[6],tempvec[3],mu);
        side_link_force(mu,sig,-Lepage,tempvec[0],tempvec[6],tempvec[1],tempvec[3]);
/* Add the P5nu vector to P3 */
        coeff = (Lepage / ThreeSt);
        for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) 
          scalar_mult_add_su3_vector((tempvec[7] + i),(tempvec[3] + i),coeff,(tempvec[7] + i));
/* Length 3 paths (Not the Naik term) */
/*Add the force F_mu the 0(2) link in the path: +     */
        if (mu <= 3) 
          u_shift_fermion(tempvec[7],tempvec[3],mu);
/* The above shift is not needed if mu is backwards */
        side_link_force(mu,sig,ThreeSt,temp_x,tempvec[7],tempvec[0],tempvec[3]);
/* Finally the OneLink and the Naik term */
/* Check if this direction is not already done */
        if (!(DirectLinks[mu] != 0)) {
          DirectLinks[mu] = 1;
/* Do only the forward terms in the Dslash */
          if (mu > 3) {
/* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
		 * shift.                                                   */
/* The one link */
            add_force_to_mom(tempvec[0],temp_x,(7 - mu),OneLink);
/* For the same reason Pnumu is the forward double link */
/* Popmu is a backward shift */
            u_shift_fermion(temp_x,tempvec[4],mu);
/* The Naik */
/* link no 1: - */
            add_force_to_mom(tempvec[1],tempvec[4],(7 - mu),-Naik);
/*Pmumumu can overwrite Popmu which is no longer needed */
            u_shift_fermion(tempvec[1],tempvec[4],(7 - mu));
/* link no 0: + */
            add_force_to_mom(tempvec[4],temp_x,(7 - mu),Naik);
          }
          else 
/* The rest of the Naik terms */
{
            u_shift_fermion(temp_x,tempvec[4],mu);
/* link no 2: + */
/* Pnumu is double backward shift */
            add_force_to_mom(tempvec[4],tempvec[1],mu,Naik);
          }
        }
/* mu */
      }
/* Here we have to do together the Naik term and the one link term */
/*sig */
  }
/* Free temporary vectors */
  free(temp_x);
  for (mu = 0; mu < 8; mu++) 
    free(tempvec[mu]);
#ifdef FFTIME
  dtime += dclock();
  if (this_node == 0) 
    printf("FFTIME:  time = %e mflops = %e\n",dtime,((((float )nflop) * volume) / ((1e6 * dtime) * (numnodes()))));
/**printf("TLENGTH: %d\n",tlength);**/
#endif
/* eo_fermion_force(version 6) */
}
#undef Pmu          
#undef Pnumu        
#undef Prhonumu     
#undef P7           
#undef P7rho        
#undef P7rhonu      
#undef P5           
#undef P3           
#undef P5nu         
#undef P3mu         
#undef Popmu        
#undef Pmumumu      
su3_matrix *forwardlink[4UL];
su3_matrix *tempmom[4UL];

void eo_fermion_force_3f(float eps,int nflav1,field_offset x1_off,int nflav2,field_offset x2_off)
{
/* note CG_solution and Dslash * solution are combined in "x_off" */
/* New version 1/21/99.  Use forward part of Dslash to get force */
/* 4/15/99 combine force from two different mass quarks, (eg 2+1flavors) */
/* see long comment at end */
/* For each link we need x_off transported from both ends of path. */
  int i;
  site *s;
  int mu;
  int nu;
  int rho;
  int sig;
  int dir;
  float coeff[2UL];
  float ferm_epsilon;
  float OneLink[2UL];
  float Lepage[2UL];
  float Naik[2UL];
  float FiveSt[2UL];
  float ThreeSt[2UL];
  float SevenSt[2UL];
  float mNaik[2UL];
  float mLepage[2UL];
  float mFiveSt[2UL];
  float mThreeSt[2UL];
  float mSevenSt[2UL];
  half_wilson_vector **Pnumu;
  half_wilson_vector **Prhonumu;
  half_wilson_vector **P7;
  half_wilson_vector **P7rho;
  half_wilson_vector **P5nu;
  half_wilson_vector **P3mu;
  half_wilson_vector **P5sig;
  half_wilson_vector **Popmu;
  half_wilson_vector **Pmumumu;
  half_wilson_vector *P3[8UL];
  half_wilson_vector *P5[8UL];
  half_wilson_vector *temp_x;
  half_wilson_vector *Pmu;
  half_wilson_vector *Pmumu;
  half_wilson_vector *temp_hw[8UL];
  msg_tag *mt[8UL];
  msg_tag *mtag[4UL];
  half_wilson_vector **hw_pt[8UL];
#ifdef FFTIME
  int nflop = 433968;
  double dtime;
  dtime = -dclock();
#endif
/* Allocate temporary vectors */
  for (mu = 0; mu < 8; mu++) {
    P3[mu] = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    if (P3[mu] == ((half_wilson_vector *)((void *)0))) {
      printf("eo_fermion_force_3f: No room for P3\n");
      terminate(1);
    }
    P5[mu] = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    if (P5[mu] == ((half_wilson_vector *)((void *)0))) {
      printf("eo_fermion_force_3f: No room for P5\n");
      terminate(1);
    }
  }
  Pmu = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
  if (Pmu == ((half_wilson_vector *)((void *)0))) {
    printf("eo_fermion_force_3f: No room for Pmu\n");
    terminate(1);
  }
  Pmumu = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
  if (Pmumu == ((half_wilson_vector *)((void *)0))) {
    printf("eo_fermion_force_3f: No room for Pmumu\n");
    terminate(1);
  }
/* Initialize message pointers */
  for (mu = 0; mu < 8; mu++) 
    mt[mu] = ((msg_tag *)((void *)0));
/* Allocate communication pointers (gen_pt type) */
  for (mu = 0; mu < 8; mu++) {
    half_wilson_vector **pt;
    pt = ((half_wilson_vector **)(malloc((sites_on_node * sizeof(half_wilson_vector *)))));
    if (pt == ((half_wilson_vector **)((void *)0))) {
      printf("eo_fermion_force_3f: No room for hw_pt\n");
      terminate(1);
    }
    hw_pt[mu] = pt;
  }
/* Double store forward gauge links */
  for (dir = 0; dir <= 3; dir++) {
    su3_matrix *pt;
    pt = ((su3_matrix *)(malloc((sites_on_node * sizeof(su3_matrix )))));
    if (pt == ((su3_matrix *)((void *)0))) {
      printf("eo_fermion_force_3f: No room for forwardlink\n");
      terminate(1);
    }
    forwardlink[dir] = pt;
  }
/* Gather forward links */
  for (dir = 0; dir <= 3; dir++) {
    mtag[dir] = start_gather(((field_offset )(((char *)(lattice[0].link + dir)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),(7 - dir),3,gen_pt[dir]);
  }
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(mtag[dir]);
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
      forwardlink[dir][i] =  *((su3_matrix *)gen_pt[dir][i]);
    }
    cleanup_gather(mtag[dir]);
  }
/* Uncompress gauge momenta */
  for (dir = 0; dir <= 3; dir++) {
    su3_matrix *pt;
    pt = ((su3_matrix *)(malloc((sites_on_node * sizeof(su3_matrix )))));
    if (pt == ((su3_matrix *)((void *)0))) {
      printf("eo_fermion_force_3f: No room for tempmom\n");
      terminate(1);
    }
    tempmom[dir] = pt;
  }
  for (dir = 0; dir <= 3; dir++) {
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
      uncompress_anti_hermitian(((s -> mom) + dir),(tempmom[dir] + i));
    }
  }
/* Path coefficients times fermion epsilon */
  ferm_epsilon = ((2.0 * (nflav1 / 4.0)) * eps);
  OneLink[0] = (act_path_coeff[0] * ferm_epsilon);
  Naik[0] = (act_path_coeff[1] * ferm_epsilon);
  mNaik[0] = -Naik[0];
  ThreeSt[0] = (act_path_coeff[2] * ferm_epsilon);
  mThreeSt[0] = -ThreeSt[0];
  FiveSt[0] = (act_path_coeff[3] * ferm_epsilon);
  mFiveSt[0] = -FiveSt[0];
  SevenSt[0] = (act_path_coeff[4] * ferm_epsilon);
  mSevenSt[0] = -SevenSt[0];
  Lepage[0] = (act_path_coeff[5] * ferm_epsilon);
  mLepage[0] = -Lepage[0];
  ferm_epsilon = ((2.0 * (nflav2 / 4.0)) * eps);
  OneLink[1] = (act_path_coeff[0] * ferm_epsilon);
  Naik[1] = (act_path_coeff[1] * ferm_epsilon);
  mNaik[1] = -Naik[1];
  ThreeSt[1] = (act_path_coeff[2] * ferm_epsilon);
  mThreeSt[1] = -ThreeSt[1];
  FiveSt[1] = (act_path_coeff[3] * ferm_epsilon);
  mFiveSt[1] = -FiveSt[1];
  SevenSt[1] = (act_path_coeff[4] * ferm_epsilon);
  mSevenSt[1] = -SevenSt[1];
  Lepage[1] = (act_path_coeff[5] * ferm_epsilon);
  mLepage[1] = -Lepage[1];
/* *************************************** */
/* Allocate temporary vectors */
  for (mu = 0; mu < 8; mu++) {
    temp_hw[mu] = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
    if (temp_hw[mu] == ((half_wilson_vector *)((void *)0))) {
      printf("eo_fermion_force_3f: No room for temp_hw\n");
      terminate(1);
    }
  }
/* copy x_off to a temporary vector */
  temp_x = ((half_wilson_vector *)(malloc((sites_on_node * sizeof(half_wilson_vector )))));
  for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
    temp_x[i].h[0] =  *((su3_vector *)(((char *)s) + x1_off));
    temp_x[i].h[1] =  *((su3_vector *)(((char *)s) + x2_off));
  }
  for (mu = 0; mu < 8; mu++) {
    u_shift_hw_fermion_np(temp_x,hw_pt[7 - mu],(7 - mu),(mt + (7 - mu)),temp_hw[7 - mu]);
/* Resolve the pointers for longer term use */
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
      Pmu[i] =  *hw_pt[7 - mu][i];
    }
    for (sig = 0; sig < 8; sig++) 
      if ((sig != mu) && (sig != (7 - mu))) {
        u_shift_hw_fermion_np(Pmu,hw_pt[sig],sig,(mt + sig),temp_hw[sig]);
/* Resolve the pointers for longer term use */
        for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
          P3[sig][i] =  *hw_pt[sig][i];
        }
        if (sig <= 3) {
/* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
          add_3f_force_to_mom_nn(P3[sig],Pmu,sig,mThreeSt);
        }
      }
    for (nu = 0; nu < 8; nu++) 
      if ((nu != mu) && (nu != (7 - mu))) {
        Pnumu = hw_pt[7 - nu];
        u_shift_hw_fermion_np(Pmu,Pnumu,(7 - nu),(mt + (7 - nu)),temp_hw[7 - nu]);
        for (sig = 0; sig < 8; sig++) 
          if ((((sig != mu) && (sig != (7 - mu))) && (sig != nu)) && (sig != (7 - nu))) {
            u_shift_hw_fermion_pp(Pnumu,hw_pt[sig],sig,(mt + sig),temp_hw[sig]);
/* Resolve the pointers for longer term use */
            for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
              P5[sig][i] =  *hw_pt[sig][i];
            }
            if (sig <= 3) {
/* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
              add_3f_force_to_mom_np(P5[sig],Pnumu,sig,FiveSt);
            }
          }
        for (rho = 0; rho < 8; rho++) 
          if ((((rho != mu) && (rho != (7 - mu))) && (rho != nu)) && (rho != (7 - nu))) {
            Prhonumu = hw_pt[7 - rho];
            u_shift_hw_fermion_pp(Pnumu,Prhonumu,(7 - rho),(mt + (7 - rho)),temp_hw[7 - rho]);
            for (sig = 0; sig < 8; sig++) 
              if ((((((sig != mu) && (sig != (7 - mu))) && (sig != nu)) && (sig != (7 - nu))) && (sig != rho)) && (sig != (7 - rho))) {
/* Length 7 paths */
                P7 = hw_pt[sig];
                u_shift_hw_fermion_pp(Prhonumu,P7,sig,(mt + sig),temp_hw[sig]);
                if (sig <= 3) {
/* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
                  add_3f_force_to_mom_pp(P7,Prhonumu,sig,mSevenSt);
                }
/* Add the force F_rho the 2(4) link in the path: +     */
                P7rho = hw_pt[rho];
                u_shift_hw_fermion_pp(P7,P7rho,rho,(mt + rho),temp_hw[rho]);
                side_link_3f_force_pppp(rho,sig,SevenSt,Pnumu,P7,Prhonumu,P7rho);
/* Add the P7rho vector to P5 */
                coeff[0] = (SevenSt[0] / FiveSt[0]);
                coeff[1] = (SevenSt[1] / FiveSt[1]);
                for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
                  scalar_mult_add_su3_vector((P5[sig][i].h + 0),(( *P7rho[i]).h + 0),coeff[0],(P5[sig][i].h + 0));
                  scalar_mult_add_su3_vector((P5[sig][i].h + 1),(( *P7rho[i]).h + 1),coeff[1],(P5[sig][i].h + 1));
                }
/* sig */
              }
/* rho */
          }
        for (sig = 0; sig < 8; sig++) 
          if ((((sig != mu) && (sig != (7 - mu))) && (sig != nu)) && (sig != (7 - nu))) {
/* Length 5 paths */
/* Add the force F_nu the 1(3) link in the path: -     */
            P5nu = hw_pt[nu];
            u_shift_hw_fermion_np(P5[sig],P5nu,nu,(mt + nu),temp_hw[nu]);
            side_link_3f_force_nnpp(nu,sig,mFiveSt,Pmu,P5[sig],Pnumu,P5nu);
/* Add the P5nu vector to P3 */
            coeff[0] = (FiveSt[0] / ThreeSt[0]);
            coeff[1] = (FiveSt[1] / ThreeSt[1]);
            for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
              scalar_mult_add_su3_vector((P3[sig][i].h + 0),(( *P5nu[i]).h + 0),coeff[0],(P3[sig][i].h + 0));
              scalar_mult_add_su3_vector((P3[sig][i].h + 1),(( *P5nu[i]).h + 1),coeff[1],(P3[sig][i].h + 1));
            }
/* sig */
          }
/* nu */
      }
/* Now the Lepage term... It is the same as 5-link paths with
	 nu=mu and FiveSt=Lepage. */
    u_shift_hw_fermion_np(Pmu,hw_pt[7 - mu],(7 - mu),(mt + (7 - mu)),temp_hw[7 - mu]);
/* Resolve the pointers for longer term use */
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
      Pmumu[i] =  *hw_pt[7 - mu][i];
    }
    for (sig = 0; sig < 8; sig++) 
      if ((sig != mu) && (sig != (7 - mu))) {
        P5sig = hw_pt[sig];
        u_shift_hw_fermion_np(Pmumu,P5sig,sig,(mt + sig),temp_hw[sig]);
        if (sig <= 3) {
/* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
          add_3f_force_to_mom_pn(P5sig,Pmumu,sig,Lepage);
        }
/* Add the force F_nu the 1(3) link in the path: -     */
        P5nu = hw_pt[mu];
        u_shift_hw_fermion_pp(P5sig,P5nu,mu,(mt + mu),temp_hw[mu]);
        side_link_3f_force_npnp(mu,sig,mLepage,Pmu,P5sig,Pmumu,P5nu);
/* Add the P5nu vector to P3 */
        coeff[0] = (Lepage[0] / ThreeSt[0]);
        coeff[1] = (Lepage[1] / ThreeSt[1]);
        for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
          scalar_mult_add_su3_vector((P3[sig][i].h + 0),(( *P5nu[i]).h + 0),coeff[0],(P3[sig][i].h + 0));
          scalar_mult_add_su3_vector((P3[sig][i].h + 1),(( *P5nu[i]).h + 1),coeff[1],(P3[sig][i].h + 1));
        }
/* Length 3 paths (Not the Naik term) */
/* Add the force F_mu the 0(2) link in the path: +     */
        if (mu <= 3) 
/* OK to clobber P5nu */
          P3mu = hw_pt[mu];
        u_shift_hw_fermion_np(P3[sig],P3mu,mu,(mt + mu),temp_hw[mu]);
/* The above shift is not needed if mu is backwards */
        side_link_3f_force_nnnp(mu,sig,ThreeSt,temp_x,P3[sig],Pmu,P3mu);
      }
/* Finally the OneLink and the Naik term */
/* Do only the forward terms in the Dslash */
    if (mu > 3) {
/* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
	   * shift.                                                   */
/* The one link */
      add_3f_force_to_mom_nn(Pmu,temp_x,(7 - mu),OneLink);
/* For the same reason Pmumu is the forward double link */
/* Popmu is a backward shift */
/* OK to clobber P3mu */
      Popmu = hw_pt[mu];
      u_shift_hw_fermion_np(temp_x,Popmu,mu,(mt + mu),temp_hw[mu]);
/* The Naik */
/* link no 1: - */
      add_3f_force_to_mom_np(Pmumu,Popmu,(7 - mu),mNaik);
/* Pmumumu can overwrite Popmu which is no longer needed */
      Pmumumu = hw_pt[7 - mu];
      u_shift_hw_fermion_np(Pmumu,Pmumumu,(7 - mu),(mt + (7 - mu)),temp_hw[7 - mu]);
/* link no 0: + */
      add_3f_force_to_mom_pn(Pmumumu,temp_x,(7 - mu),Naik);
    }
    else 
/* The rest of the Naik terms */
{
/* OK to clobber P3mu */
      Popmu = hw_pt[mu];
      u_shift_hw_fermion_np(temp_x,Popmu,mu,(mt + mu),temp_hw[mu]);
/* link no 2: + */
/* Pmumu is double backward shift */
      add_3f_force_to_mom_pn(Popmu,Pmumu,mu,Naik);
    }
/* Here we have to do together the Naik term and the one link term */
/* mu */
  }
/* Repack momenta */
  for (dir = 0; dir <= 3; dir++) {
    for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
      make_anti_hermitian((tempmom[dir] + i),((s -> mom) + dir));
    }
  }
/* Free temporary vectors */
  free(temp_x);
  free(Pmu);
  free(Pmumu);
  for (mu = 0; mu < 8; mu++) {
    free(P3[mu]);
    free(P5[mu]);
    free(hw_pt[mu]);
    free(temp_hw[mu]);
  }
  for (dir = 0; dir <= 3; dir++) {
    free(forwardlink[dir]);
    free(tempmom[dir]);
  }
/* Cleanup gathers */
  for (mu = 0; mu < 8; mu++) 
    if (mt[mu] != ((msg_tag *)((void *)0))) 
      cleanup_gather(mt[mu]);
#ifdef FFTIME
  dtime += dclock();
  if (this_node == 0) 
    printf("FFTIME:  time = %e mflops = %e\n",dtime,((((float )nflop) * volume) / ((1e6 * dtime) * (numnodes()))));
/**printf("TLENGTH: %d\n",tlength);**/
#endif
/* eo_fermion_force_3f */
}
#undef Pmu          
#undef Pnumu        
#undef Prhonumu     
#undef P7           
#undef P7rho        
#undef P7rhonu      
#undef P5           
#undef P3           
#undef P5nu         
#undef P3mu         
#undef Popmu        
#undef Pmumumu      
#endif /* ASQ_OPTIMIZED_FORCE */
#ifdef  ASQ_OPTIMIZED_FATTENING   /* Asqtad action only, "_fn" executables */
#ifndef FN
#endif
#ifdef DM_DU0
#else

void compute_gen_staple(field_offset staple,int mu,int nu,field_offset link,float coef)
{
#endif
  msg_tag *mtag0;
  msg_tag *mtag1;
  su3_matrix *tempmat;
  register site *s;
  register int i;
/* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any su3_matrix. The result is saved in staple.
    if staple==NULL_FP then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */
/* Upper staple */
  mtag0 = start_gather(link,(sizeof(su3_matrix )),nu,3,gen_pt[0]);
  mtag1 = start_gather(((field_offset )(((char *)(lattice[0].link + nu)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),mu,3,gen_pt[1]);
  wait_gather(mtag0);
  wait_gather(mtag1);
/* Save the staple */
  if (staple != -1) {
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      su3_matrix tmat1;
// mult_su3_na( (su3_matrix *)gen_pt[0][i],
      asr_su3_na(((su3_matrix *)gen_pt[0][i]),((su3_matrix *)gen_pt[1][i]),&tmat1);
// mult_su3_nn( &(s->link[nu]), &tmat1, (su3_matrix *)F_PT(s,staple) );
      asr_su3_nn(((s -> link) + nu),&tmat1,((su3_matrix *)(((char *)s) + staple)));
      s++;
    }
  }
  else 
/* No need to save the staple. Add it to the fatlinks */
{
#if 1
    su3_matrix *pTMAT1;
    su3_matrix *pTMAT2;
    su3_matrix **pFAT1;
    pTMAT1 = ((su3_matrix *)(malloc((sizeof(su3_matrix ) * sites_on_node))));
    pTMAT2 = ((su3_matrix *)(malloc((sizeof(su3_matrix ) * sites_on_node))));
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
// mult_su3_na( (su3_matrix *)gen_pt[0][i],
      asr_su3_na(((su3_matrix *)gen_pt[0][i]),((su3_matrix *)gen_pt[1][i]),(pTMAT1 + i));
// mult_su3_nn( &(s->link[nu]), &(pTMAT1[i]), &(pTMAT2[i]) );
      asr_su3_nn(((s -> link) + nu),(pTMAT1 + i),(pTMAT2 + i));
      s++;
    }
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
#ifdef DSLASH_TMP_LINKS
      scalar_mult_add_su3_matrix((t_fatlink + ((4 * i) + mu)),(pTMAT2 + i),coef,(t_fatlink + ((4 * i) + mu)));
#else
#endif
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
#endif
      s++;
    }
    free(pTMAT1);
    free(pTMAT2);
#else
#ifdef DSLASH_TMP_LINKS
#else
#endif
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
#endif
#endif
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);
/* lower staple */
  posix_memalign((&tempmat),256,(sites_on_node * sizeof(su3_matrix )));
  mtag0 = start_gather(((field_offset )(((char *)(lattice[0].link + nu)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),mu,3,gen_pt[0]);
  wait_gather(mtag0);
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    su3_matrix tmat1;
    asr_su3_an(((s -> link) + nu),((su3_matrix *)(((char *)s) + link)),&tmat1);
    asr_su3_nn(&tmat1,((su3_matrix *)gen_pt[0][i]),(tempmat + i));
    s++;
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_from_temp(tempmat,(sizeof(su3_matrix )),(7 - nu),3,gen_pt[0]);
  wait_gather(mtag0);
/* Save the staple */
  if (staple != -1) {
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      register su3_matrix *fat1;
      add_su3_matrix(((su3_matrix *)(((char *)s) + staple)),((su3_matrix *)gen_pt[0][i]),((su3_matrix *)(((char *)s) + staple)));
#ifdef DSLASH_TMP_LINKS
      fat1 = (t_fatlink + ((4 * i) + mu));
#else
#endif
      scalar_mult_add_su3_matrix(fat1,((su3_matrix *)(((char *)s) + staple)),coef,fat1);
#ifdef DM_DU0
#ifdef DSLASH_TMP_LINKS
#else
#endif
#endif
      s++;
    }
  }
  else 
/* No need to save the staple. Add it to the fatlinks */
{
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      register su3_matrix *fat1;
#ifdef DSLASH_TMP_LINKS
      fat1 = (t_fatlink + ((4 * i) + mu));
#else
#endif
      scalar_mult_add_su3_matrix(fat1,((su3_matrix *)gen_pt[0][i]),coef,fat1);
      s++;
    }
  }
  free(tempmat);
  cleanup_gather(mtag0);
}
#endif  /* ASQ_OPTIMIZED_FATTENING   */
#ifdef  ASQ_OPTIMIZED_FORCE
/*   Covariant shift of the src fermion field in the direction dir  *
 *  by one unit. The result is stored in dest.                       */

void u_shift_fermion(su3_vector *src,su3_vector *dest,int dir)
{
  su3_vector *tmpvec;
  msg_tag *mtag;
  register site *s;
  register int i;
/* forward shift */
  if (dir <= 3) {
    mtag = start_gather_from_temp(src,(sizeof(su3_vector )),dir,3,gen_pt[0]);
    wait_gather(mtag);
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_su3_mat_vec(((s -> link) + dir),((su3_vector *)gen_pt[0][i]),(dest + i));
      s++;
    }
    cleanup_gather(mtag);
  }
  else 
/* backward shift */
{
    tmpvec = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_adj_su3_mat_vec(((s -> link) + (7 - dir)),(src + i),(tmpvec + i));
      s++;
    }
    mtag = start_gather_from_temp(tmpvec,(sizeof(su3_vector )),dir,3,gen_pt[0]);
    wait_gather(mtag);
/* copy the gen_pt to the dest */
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      dest[i] =  *((su3_vector *)gen_pt[0][i]);
      s++;
    }
    cleanup_gather(mtag);
    free(tmpvec);
  }
}
/*  Covariant shift of the src half wilson fermion field in the  *
 * direction dir   by one unit.  The result is stored in dest pointer.    */

void u_shift_hw_fermion_np(half_wilson_vector *src,half_wilson_vector **dest_pt,int dir,msg_tag **mtag,half_wilson_vector *tmpvec)
{
  site *s;
  int i;
/* forward shift */
  if (dir <= 3) {
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_su3_mat_hwvec((forwardlink[dir] + i),(src + i),(tmpvec + i));
      s++;
    }
  }
  else 
/* backward shift */
{
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_adj_su3_mat_hwvec(((s -> link) + (7 - dir)),(src + i),(tmpvec + i));
      s++;
    }
  }
  if ( *mtag == ((msg_tag *)((void *)0))) 
     *mtag = start_gather_from_temp(tmpvec,(sizeof(half_wilson_vector )),dir,3,((char **)dest_pt));
  else 
    restart_gather_from_temp(tmpvec,(sizeof(half_wilson_vector )),dir,3,((char **)dest_pt), *mtag);
  wait_gather( *mtag);
}
/*  Covariant shift of the src half wilson fermion field in the  *
 * direction dir   by one unit. The result is stored in dest.   
 * This version shifts from a list of pointers */

void u_shift_hw_fermion_pp(half_wilson_vector **src_pt,half_wilson_vector **dest_pt,int dir,msg_tag **mtag,half_wilson_vector *tmpvec)
{
  site *s;
  int i;
/* forward shift */
  if (dir <= 3) {
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_su3_mat_hwvec((forwardlink[dir] + i),src_pt[i],(tmpvec + i));
      s++;
    }
  }
  else 
/* backward shift */
{
    s = lattice;
// #pragma simd
    for (i = 0; i < sites_on_node; i++) {
      mult_adj_su3_mat_hwvec(((s -> link) + (7 - dir)),src_pt[i],(tmpvec + i));
      s++;
    }
  }
  if ( *mtag == ((msg_tag *)((void *)0))) 
     *mtag = start_gather_from_temp(tmpvec,(sizeof(half_wilson_vector )),dir,3,((char **)dest_pt));
  else 
    restart_gather_from_temp(tmpvec,(sizeof(half_wilson_vector )),dir,3,((char **)dest_pt), *mtag);
  wait_gather( *mtag);
}
/* Add in contribution to the force */
/* Put antihermitian traceless part into momentum */

void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,float coeff)
{
  register site *s;
  register int i;
  register float tmp_coeff;
  su3_matrix tmat;
  su3_matrix tmat2;
  if (dir > 3) {
    dir = (7 - dir);
    coeff = -coeff;
  }
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    if ((s -> parity) == 1) 
      tmp_coeff = -coeff;
    else 
      tmp_coeff = coeff;
    uncompress_anti_hermitian(((s -> mom) + dir),&tmat2);
    su3_projector((back + i),(forw + i),&tmat);
    scalar_mult_add_su3_matrix(&tmat2,&tmat,tmp_coeff,&tmat2);
    make_anti_hermitian(&tmat2,((s -> mom) + dir));
    s++;
  }
}
/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */

void add_3f_force_to_mom_nn(half_wilson_vector *back,half_wilson_vector *forw,int dir,float coeff[2UL])
{
  register site *s;
  register int i;
  float tmp_coeff[2UL];
  su3_matrix tmat;
  su3_matrix *tmat2;
  if (dir > 3) {
    dir = (7 - dir);
    coeff[0] = -coeff[0];
    coeff[1] = -coeff[1];
  }
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    if ((s -> parity) == 1) {
      tmp_coeff[0] = -coeff[0];
      tmp_coeff[1] = -coeff[1];
    }
    else {
      tmp_coeff[0] = coeff[0];
      tmp_coeff[1] = coeff[1];
    }
    tmat2 = (tempmom[dir] + i);
    su3_projector((back[i].h + 0),(forw[i].h + 0),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[0],tmat2);
    su3_projector((back[i].h + 1),(forw[i].h + 1),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[1],tmat2);
    s++;
  }
}
/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
/* This variant takes a list of pointers for back and forw */

void add_3f_force_to_mom_np(half_wilson_vector *back,half_wilson_vector **forw_pt,int dir,float (coeff)[2UL])
{
  register site *s;
  register int i;
  float tmp_coeff[2UL];
  su3_matrix tmat;
  su3_matrix *tmat2;
  if (dir > 3) {
    dir = (7 - dir);
    coeff[0] = -coeff[0];
    coeff[1] = -coeff[1];
  }
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    if ((s -> parity) == 1) {
      tmp_coeff[0] = -coeff[0];
      tmp_coeff[1] = -coeff[1];
    }
    else {
      tmp_coeff[0] = coeff[0];
      tmp_coeff[1] = coeff[1];
    }
    tmat2 = (tempmom[dir] + i);
    su3_projector((back[i].h + 0),(( *forw_pt[i]).h + 0),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[0],tmat2);
    su3_projector((back[i].h + 1),(( *forw_pt[i]).h + 1),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[1],tmat2);
    s++;
  }
}
/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
/* This variant takes a list of pointers for back and forw */

void add_3f_force_to_mom_pn(half_wilson_vector **back_pt,half_wilson_vector *forw,int dir,float (coeff)[2UL])
{
  register site *s;
  register int i;
  float tmp_coeff[2UL];
  su3_matrix tmat;
  su3_matrix *tmat2;
  if (dir > 3) {
    dir = (7 - dir);
    coeff[0] = -coeff[0];
    coeff[1] = -coeff[1];
  }
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    if ((s -> parity) == 1) {
      tmp_coeff[0] = -coeff[0];
      tmp_coeff[1] = -coeff[1];
    }
    else {
      tmp_coeff[0] = coeff[0];
      tmp_coeff[1] = coeff[1];
    }
    tmat2 = (tempmom[dir] + i);
    su3_projector((( *back_pt[i]).h + 0),(forw[i].h + 0),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[0],tmat2);
    su3_projector((( *back_pt[i]).h + 1),(forw[i].h + 1),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[1],tmat2);
    s++;
  }
}
/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
/* This variant takes a list of pointers for back and forw */

void add_3f_force_to_mom_pp(half_wilson_vector **back_pt,half_wilson_vector **forw_pt,int dir,float (coeff)[2UL])
{
  register site *s;
  register int i;
  float tmp_coeff[2UL];
  su3_matrix tmat;
  su3_matrix *tmat2;
  if (dir > 3) {
    dir = (7 - dir);
    coeff[0] = -coeff[0];
    coeff[1] = -coeff[1];
  }
  s = lattice;
// #pragma simd
  for (i = 0; i < sites_on_node; i++) {
    if ((s -> parity) == 1) {
      tmp_coeff[0] = -coeff[0];
      tmp_coeff[1] = -coeff[1];
    }
    else {
      tmp_coeff[0] = coeff[0];
      tmp_coeff[1] = coeff[1];
    }
    tmat2 = (tempmom[dir] + i);
    su3_projector((( *back_pt[i]).h + 0),(( *forw_pt[i]).h + 0),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[0],tmat2);
    su3_projector((( *back_pt[i]).h + 1),(( *forw_pt[i]).h + 1),&tmat);
    scalar_mult_add_su3_matrix(tmat2,&tmat,tmp_coeff[1],tmat2);
    s++;
  }
}
/*  This routine is needed in order to add the force on the side link *
 * of the paths in the Asq and Asqtad actions. It gets as inputs the  *
 * direction mu of the side link and the direction nu of the Dslash   *
 * term we are dealing with. Then it takes also 4 fermion fields:     *
 * Path: the piece of the path with no hop in the nu or mu direction  *
 * Path_nu: the piece of the path with a hop in nu  but not in mu     *
 * Path_mu: is Path times the link mu                                 *
 * Path_numu: is Path_nu times the link mu                            */

void side_link_force(int mu,int nu,float coeff,su3_vector *Path,su3_vector *Path_nu,su3_vector *Path_mu,su3_vector *Path_numu)
{
  if (mu <= 3) {
/*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
    if (nu <= 3) 
      add_force_to_mom(Path_numu,Path,mu,coeff);
    else 
/* ? extra - */
      add_force_to_mom(Path,Path_numu,(7 - mu),-coeff);
  }
  else 
/*GOES_BACKWARDS(mu)*/
{
/* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */
    if (nu <= 3) 
/* ? extra - */
      add_force_to_mom(Path_nu,Path_mu,mu,-coeff);
    else 
      add_force_to_mom(Path_mu,Path_nu,(7 - mu),coeff);
  }
}
/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */

void side_link_3f_force_nnpp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector *Path_nu,half_wilson_vector **Path_mu,half_wilson_vector **Path_numu)
{
  float m_coeff[2UL];
  m_coeff[0] = -coeff[0];
  m_coeff[1] = -coeff[1];
  if (mu <= 3) {
/*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
    if (nu <= 3) 
      add_3f_force_to_mom_pn(Path_numu,Path,mu,coeff);
    else 
/* ? extra - */
      add_3f_force_to_mom_np(Path,Path_numu,(7 - mu),m_coeff);
  }
  else 
/*GOES_BACKWARDS(mu)*/
{
/* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */
    if (nu <= 3) 
/* ? extra - */
      add_3f_force_to_mom_np(Path_nu,Path_mu,mu,m_coeff);
    else 
      add_3f_force_to_mom_pn(Path_mu,Path_nu,(7 - mu),coeff);
  }
}
/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */

void side_link_3f_force_nnnp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector *Path_nu,half_wilson_vector *Path_mu,half_wilson_vector **Path_numu)
{
  float m_coeff[2UL];
  m_coeff[0] = -coeff[0];
  m_coeff[1] = -coeff[1];
  if (mu <= 3) {
/*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
    if (nu <= 3) 
      add_3f_force_to_mom_pn(Path_numu,Path,mu,coeff);
    else 
/* ? extra - */
      add_3f_force_to_mom_np(Path,Path_numu,(7 - mu),m_coeff);
  }
  else 
/*GOES_BACKWARDS(mu)*/
{
/* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */
    if (nu <= 3) 
/* ? extra - */
      add_3f_force_to_mom_nn(Path_nu,Path_mu,mu,m_coeff);
    else 
      add_3f_force_to_mom_nn(Path_mu,Path_nu,(7 - mu),coeff);
  }
}
/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */

void side_link_3f_force_npnp(int mu,int nu,float coeff[2UL],half_wilson_vector *Path,half_wilson_vector **Path_nu,half_wilson_vector *Path_mu,half_wilson_vector **Path_numu)
{
  float m_coeff[2UL];
  m_coeff[0] = -coeff[0];
  m_coeff[1] = -coeff[1];
  if (mu <= 3) {
/*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
    if (nu <= 3) 
      add_3f_force_to_mom_pn(Path_numu,Path,mu,coeff);
    else 
/* ? extra - */
      add_3f_force_to_mom_np(Path,Path_numu,(7 - mu),m_coeff);
  }
  else 
/*GOES_BACKWARDS(mu)*/
{
/* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */
    if (nu <= 3) 
/* ? extra - */
      add_3f_force_to_mom_pn(Path_nu,Path_mu,mu,m_coeff);
    else 
      add_3f_force_to_mom_np(Path_mu,Path_nu,(7 - mu),coeff);
  }
}
/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */

void side_link_3f_force_pppp(int mu,int nu,float coeff[2UL],half_wilson_vector **Path,half_wilson_vector **Path_nu,half_wilson_vector **Path_mu,half_wilson_vector **Path_numu)
{
  float m_coeff[2UL];
  m_coeff[0] = -coeff[0];
  m_coeff[1] = -coeff[1];
  if (mu <= 3) {
/*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
    if (nu <= 3) 
      add_3f_force_to_mom_pp(Path_numu,Path,mu,coeff);
    else 
/* ? extra - */
      add_3f_force_to_mom_pp(Path,Path_numu,(7 - mu),m_coeff);
  }
  else 
/*GOES_BACKWARDS(mu)*/
{
/* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */
    if (nu <= 3) 
/* ? extra - */
      add_3f_force_to_mom_pp(Path_nu,Path_mu,mu,m_coeff);
    else 
      add_3f_force_to_mom_pp(Path_mu,Path_nu,(7 - mu),coeff);
  }
}
#endif  /* ASQ_OPTIMIZED_FORCE   */
/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at
< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.
The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X
 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.
Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:
 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 
	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.
The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:
	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)
	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)
	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)
As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.
Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.
Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.
*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.
Examples:
  Fat Link:
    Positive:	X-------X
    Negative     --------
	 	|	|
		|	|
		X	X
  Naik connection, smeared
    Positive:	X-------x-------x-------X
    Negative:	---------
		|	|
		|	|
		X	x-------x-------X
    Positive:	--------x--------
		|		|
		|		|
		X		x-------X
    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/
/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y
		X
		|
		|
		X
		|
		|
	X---->--X
  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
