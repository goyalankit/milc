/********************** path_product.c ***************************/
/* MIMD version 6 */
/* Compute product of links along a specified path */
/* On return lattice[i].tempmat1 contains the product for the path
   ENDING at site i.  e.g., for the 1-link path XUP,
   lattice[i].tempmat1 = lattice[i-\hat x][XUP] */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/prefetch.h"
#define FETCH_UP 1
/* LOOPEND is required now -CD */
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

void path_product(const int *dir,const int length)
{
  register int i;
  register site *s;
  msg_tag *mtag0;
  su3_matrix *tempmat2t;
  su3_matrix *tempmat3t;
  int j;
/* a forward step leaves the answer in gen_pt[0], which points into
	link, tempmat1 or tempmat2, and backwards step in tempmat1 or tempmat2,
	After a forwards step, need to wait and clean a gather.
	  STEP	leaves answer in
	  even # forward	gen_pt[0]->tempmat1  (gen_pt[0]->link for step 0
	  even # backward	tempmat1
	  odd  # forward	gen_pt[0]->tempmat2
	  odd  # backward	tempmat2
	*/
/* Trivial path case */
  if (length == 0) {{
      register int loopend;
      loopend = sites_on_node;
      for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
        clear_su3mat(&s -> tempmat1);
        s -> tempmat1.e[0][0].real = (s -> tempmat1.e[1][1].real = (s -> tempmat1.e[2][2].real = 1.));
      }
    }
    return ;
  }
/* allocate temporary space */
  posix_memalign((&tempmat3t),256,(sites_on_node * sizeof(su3_matrix )));
  posix_memalign((&tempmat2t),256,(sites_on_node * sizeof(su3_matrix )));
/* j=0 */
  if (dir[0] <= 3) {
    mtag0 = start_gather(((field_offset )(((char *)(lattice[0].link + dir[0])) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),(7 - dir[0]),3,gen_pt[0]);
  }
  else 
/* if GOES_BACKWARDS(dir[0]) */
{{
      register int loopend;
      loopend = sites_on_node;
      for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
        if (i < (loopend - 1)) {;
        }
        su3_adjoint(((s -> link) + (7 - dir[0])),&s -> tempmat1);
      }
    }
  }
  for (j = 1; j < length; j++) {
    if ((j % 2) == 1) {
      if (dir[j] <= 3) {
        if (dir[j - 1] <= 3) {
          wait_gather(mtag0);
{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;;
              }
              mult_su3_nn(((su3_matrix *)gen_pt[0][i]),((s -> link) + dir[j]),(tempmat2t + i));
            }
          }
          cleanup_gather(mtag0);
        }
        else 
/* last link was backwards */
{{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;;
              }
              mult_su3_nn(&s -> tempmat1,((s -> link) + dir[j]),(tempmat2t + i));
            }
          }
        }
        mtag0 = start_gather_from_temp(tempmat2t,(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
/* for GOES_FORWARDS */
      }
      else 
/* GOES_BACKWARDS(dir[j]), which is an odd numbered step */
{
        if (dir[j - 1] <= 3) {
          wait_gather(mtag0);
{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;
              }
              su3mat_copy(((su3_matrix *)gen_pt[0][i]),(tempmat3t + i));
            }
          }
          cleanup_gather(mtag0);
          mtag0 = start_gather_from_temp(tempmat3t,(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
        }
        else 
/*last step was backwards */
{
          mtag0 = start_gather(((field_offset )(((char *)(&lattice[0].tempmat1)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
        }
        wait_gather(mtag0);
{
          register int loopend;
          loopend = sites_on_node;
          for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
            if (i < (loopend - 1)) {;;
            }
            mult_su3_na(((su3_matrix *)gen_pt[0][i]),((s -> link) + (7 - dir[j])),(tempmat2t + i));
          }
        }
        cleanup_gather(mtag0);
/* end for GOES_BACKWARDS */
      }
/* end for j=odd */
    }
    else 
/* j=even */
{
      if (dir[j] <= 3) {
        if (dir[j - 1] <= 3) {
          wait_gather(mtag0);
{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;;
              }
              mult_su3_nn(((su3_matrix *)gen_pt[0][i]),((s -> link) + dir[j]),&s -> tempmat1);
            }
          }
          cleanup_gather(mtag0);
        }
        else 
/* last link goes backwards */
{{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;;
              }
              mult_su3_nn((tempmat2t + i),((s -> link) + dir[j]),&s -> tempmat1);
            }
          }
        }
        mtag0 = start_gather(((field_offset )(((char *)(&lattice[0].tempmat1)) - ((char *)(lattice + 0)))),(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
/* for GOES_FORWARDS */
      }
      else 
/* GOES_BACKWARDS(dir[j]) */
{
        if (dir[j - 1] <= 3) {
          wait_gather(mtag0);
{
            register int loopend;
            loopend = sites_on_node;
            for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
              if (i < (loopend - 1)) {;;
              }
              su3mat_copy(((su3_matrix *)gen_pt[0][i]),(tempmat3t + i));
            }
          }
          cleanup_gather(mtag0);
          mtag0 = start_gather_from_temp(tempmat3t,(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
        }
        else 
/* last step was backwards */
{
          mtag0 = start_gather_from_temp(tempmat2t,(sizeof(su3_matrix )),(7 - dir[j]),3,gen_pt[0]);
        }
        wait_gather(mtag0);
{
          register int loopend;
          loopend = sites_on_node;
          for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
            if (i < (loopend - 1)) {;;
            }
            mult_su3_na(((su3_matrix *)gen_pt[0][i]),((s -> link) + (7 - dir[j])),&s -> tempmat1);
          }
        }
        cleanup_gather(mtag0);
/* for GOES_BACKWARDS */
      }
/* for j=even */
    }
/* j=link in loop */
  }
/* Want to end in tempmat1 */
/* last step was odd */
  if ((length % 2) == 0) {
    if (dir[length - 1] <= 3) {
      wait_gather(mtag0);
{
        register int loopend;
        loopend = sites_on_node;
        for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
          if (i < (loopend - 1)) {;
          }
          su3mat_copy(((su3_matrix *)gen_pt[0][i]),&s -> tempmat1);
        }
      }
      cleanup_gather(mtag0);
    }
    else {{
        register int loopend;
        loopend = sites_on_node;
        for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
          if (i < (loopend - 1)) {;
          }
          su3mat_copy((tempmat2t + i),&s -> tempmat1);
        }
      }
    }
  }
  else 
/* odd length path */
{
    if (dir[length - 1] <= 3) {
      wait_gather(mtag0);
{
        register int loopend;
        loopend = sites_on_node;
        for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
          if (i < (loopend - 1)) {;
          }
          su3mat_copy(((su3_matrix *)gen_pt[0][i]),(tempmat3t + i));
        }
      }
      cleanup_gather(mtag0);
{
        register int loopend;
        loopend = sites_on_node;
        for (((i = 0) , (s = lattice)); i < loopend; (i++ , s++)) {
          if (i < (loopend - 1)) {;
          }
          su3mat_copy((tempmat3t + i),&s -> tempmat1);
        }
      }
    }
    else {
    }
  }
  free(tempmat3t);
  free(tempmat2t);
/* path */
}
#ifdef N_SUBL32
/* code from symanzik_sl32/dsdu_qhb.c ****************************/
/* U.M. Heller August 1997 */
/* This is a modification of "path_product" from gauge_stuff.c
   which works only on one sublattice. */
/* A forward step leaves the answer in gen_pt[0], which points into
       link, tempmat1 or tempmat2, and backwards step in tempmat1 or tempmat2.
       After a forwards step, need to wait and clean a gather.
	STEP			leaves answer in
	even # forward		gen_pt[0]->tempmat1 (gen_pt[0]->link for step 0)
	even # backward		tempmat1
	odd  # forward		gen_pt[0]->tempmat2
	odd  # backward		tempmat2
    */
/* allocate temporary space */
/* j=0 */
/* if GOES_BACKWARDS(dir[0]) */
/* last link was backwards */
/* for GOES_FORWARDS */
/* GOES_BACKWARDS(dir[j]), which is an odd numbered step */
/*last step was backwards */
/* end for GOES_BACKWARDS */
/* end for j=odd */
/* j=even */
/* last link was backwards */
/* for GOES_FORWARDS */
/* GOES_BACKWARDS(dir[j]), which is an even numbered step */
/*last step was backwards */
/* end for GOES_BACKWARDS */
/* end for j=even */
/* j=link in loop */
/* Want to end in tempmat1 */
/* last step was odd */
/* odd length path: last step was even */
/* path_prod_subl */
#endif /* N_SUBL32 */
