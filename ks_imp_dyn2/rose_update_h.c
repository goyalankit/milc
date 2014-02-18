/****** update_h.c  -- ******************/
/* updates momentum matrices for improved action */
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*
* MIMD version 6 */
#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h(float eps)
{
/* gauge field force */
  rephase(0);
  imp_gauge_force(eps,((field_offset )(((char *)(&lattice[0].mom)) - ((char *)(lattice + 0)))));
  rephase(1);
/* fermionic force */
/* First compute M*xxx in temporary vector xxx_odd */
/* See long comment at end of file */
/* The diagonal term in M doesn't matter */
  dslash_fn(((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),1);
  dslash_fn(((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))),((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))),1);
/**
    eo_fermion_force( eps, nflavors1, F_OFFSET(xxx1) );
    eo_fermion_force( eps, nflavors2, F_OFFSET(xxx2) );
**/
/**/
  eo_fermion_force_3f(eps,nflavors1,((field_offset )(((char *)(&lattice[0].xxx1)) - ((char *)(lattice + 0)))),nflavors2,((field_offset )(((char *)(&lattice[0].xxx2)) - ((char *)(lattice + 0)))));
/**/
/* update_h */
}
