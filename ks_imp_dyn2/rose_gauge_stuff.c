/****** gauge_stuff.c  -- ******************/
/* MIMD version 6 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file */
/**#define GFTIME**/
/* For timing gauge force calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#ifdef LOOPEND
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }
#else
#define END_LOOP        /* define it to be nothing */
#endif
#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
void printpath(int *path,int length);
#define GAUGE_ACTION_PART1
/* defines NREPS NLOOP MAX_LENGTH MAX_NUM */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART1
char gauge_action_description[128UL];
int gauge_action_nloops = 3;
int gauge_action_nreps = 1;
/* lengths of various kinds of loops */
int loop_length[3UL];
/* number of rotations/reflections  for each kind */
int loop_num[3UL];
/* table of directions, 1 for each kind of loop */
int loop_ind[3UL][6UL];
/* table of directions, for each rotation and reflection of each kind of
	loop.  tabulated with "canonical" starting point and direction. */
int loop_table[3UL][16UL][6UL];
/* table of coefficients in action, for various "representations" (actually,
	powers of the trace) */
float loop_coeff[3UL][1UL];
/* for each rotation/reflection, an integer distinct for each starting
	point, or each cyclic permutation of the links */
int loop_char[16UL];
/* for each kind of loop for each rotation/reflection, the expectation
	value of the loop */
double loop_expect[3UL][1UL][16UL];
/* Make table of loops in action */

void make_loop_table()
{
  int perm[8UL];
  int pp[8UL];
  int ir[4UL];
  int length;
  int iloop;
  int i;
  int j;
  int chr;
  int vec[6UL];
  int count;
  int flag;
  void char_num(int *dig,int *chr,int length);
#define GAUGE_ACTION_PART2
/* defines all loops and their coefficients */
#include <gauge_action.h>
#undef GAUGE_ACTION_PART2
  for (iloop = 0; iloop < 3; iloop++) {
    length = loop_length[iloop];
    count = 0;
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
                        vec[j] = pp[loop_ind[iloop][j]];
                      char_num(vec,&chr,length);
                      flag = 0;
/* check if it's a new set: */
                      for (j = 0; j < count; j++) 
                        if (chr == loop_char[j]) 
                          flag = 1;
                      if (flag == 0) {
                        loop_char[count] = chr;
                        for (j = 0; j < length; j++) 
                          loop_table[iloop][count][j] = vec[j];
                        count++;
/**node0_printf("ADD LOOP: "); printpath( vec, length );**/
                      }
                      if (count > 16) {
                        if (this_node == 0) 
                          printf("OOPS: MAX_NUM too small\n");
                        exit(0);
                      }
                      loop_num[iloop] = count;
/* end reflection*/
                    }
/* end permutation if block */
            }
/* end permutation */
          }
/* end iloop */
  }
/* print out the loop coefficients */
  if (this_node == 0) 
    printf("loop coefficients: nloop rep loop_coeff  multiplicity\n");
  for (i = 0; i < 1; i++) 
    for (j = 0; j < 3; j++) {
      if (this_node == 0) 
        printf("                    %d %d      %e     %d\n",j,i,loop_coeff[j][i],loop_num[j]);
    }
/* make_loop_table */
}
/* find a number uniquely identifying the cyclic permutation of a path,
   or the starting point on the path.  Backwards paths are considered
   equivalent here, so scan those too. */

void char_num(int *dig,int *chr,int length)
{
  int j;
  int bdig[6UL];
  int tenl;
  int newv;
  int old;
/* "dig" is array of directions.  "bdig" is array of directions for
	backwards path. */
  tenl = 1;
  for (j = 0; j < (length - 1); j++) 
    tenl = (tenl * 10);
   *chr = dig[length - 1];
  for (j = (length - 2); j >= 0; j--) 
     *chr = (( *chr * 10) + dig[j]);
/* forward*/
  old =  *chr;
  for (j = (length - 1); j >= 1; j--) {
    newv = (old - (tenl * dig[j]));
    newv = ((newv * 10) + dig[j]);
    if (newv <  *chr) 
       *chr = newv;
    old = newv;
  }
/* backward*/
  for (j = 0; j < length; j++) 
    bdig[j] = (7 - dig[(length - j) - 1]);
  old = bdig[length - 1];
  for (j = (length - 2); j >= 0; j--) 
    old = ((old * 10) + bdig[j]);
  if (old <  *chr) 
     *chr = old;
  for (j = (length - 1); j >= 1; j--) {
    newv = (old - (tenl * bdig[j]));
    newv = ((newv * 10) + bdig[j]);
    if (newv <  *chr) 
       *chr = newv;
    old = newv;
  }
/* char_num */
}

double imp_gauge_action()
{
  register int i;
  int rep;
  register site *s;
  complex trace;
  double g_action;
  double action;
  double act2;
  double total_action;
  int length;
/* these are for loop_table  */
  int ln;
  int iloop;
  g_action = 0.0;
/* gauge action */
  for (iloop = 0; iloop < 3; iloop++) {
    length = loop_length[iloop];
/* loop over rotations and reflections */
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      path_product(loop_table[iloop][ln],length);
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        trace = trace_su3(&s -> tempmat1);
        action = (3.0 - ((double )trace.real));
/* need the "3 -" for higher characters */
        total_action = (((double )loop_coeff[iloop][0]) * action);
        act2 = action;
        for (rep = 1; rep < 1; rep++) {
          act2 *= action;
          total_action += (((double )loop_coeff[iloop][rep]) * act2);
        }
        g_action += total_action;
/* sites */
      }
/* ln */
    }
/* iloop */
  }
  g_doublesum(&g_action);
  return g_action;
/* imp_gauge_action */
}
/* update the momenta with the gauge force */

void imp_gauge_force(float eps,field_offset mom_off)
{
  register int i;
  register int dir;
  register site *st;
  su3_matrix tmat1;
  su3_matrix tmat2;
  register float eb3;
  register anti_hermitmat *momentum;
/* For Symanzik1 action */
  int nflop = 153004;
#ifdef GFTIME
  double dtime;
#endif
  int j;
  int k;
  int dirs[6UL];
  int length;
  int path_dir[6UL];
  int path_length;
  int ln;
  int iloop;
  float action;
  float act2;
  float new_term;
  int ncount;
#ifdef GFTIME
  dtime = -dclock();
#endif
  eb3 = ((eps * beta) / 3.0);
/* Loop over directions, update mom[dir] */
  for (dir = 0; dir <= 3; dir++) {
    for (((i = 0) , (st = lattice)); i < sites_on_node; (i++ , st++)) 
      for (j = 0; j < 3; j++) 
        for (k = 0; k < 3; k++) {
          st -> staple.e[j][k] = cmplx(0.0,0.0);
        }
    ncount = 0;
    for (iloop = 0; iloop < 3; iloop++) {
      length = loop_length[iloop];
      for (ln = 0; ln < loop_num[iloop]; ln++) {
/**printf("UPD:  "); printpath( loop_table[iloop][ln], length );**/
/* set up dirs.  we are looking at loop starting in "XUP"
		   direction, rotate so it starts in "dir" direction. */
        for (k = 0; k < length; k++) {
          if (loop_table[iloop][ln][k] <= 3) {
            dirs[k] = ((dir + loop_table[iloop][ln][k]) % 4);
          }
          else {
            dirs[k] = (7 - ((dir + (7 - loop_table[iloop][ln][k])) % 4));
          }
        }
/* generalized "staple" */
        path_length = (length - 1);
/* check for links in direction of momentum to be
		   updated, each such link gives a contribution. Note
		   the direction of the path - opposite the link. */
        for (k = 0; k < length; k++) 
          if ((dirs[k] == dir) || (dirs[k] == (7 - dir))) {
            if (dirs[k] <= 3) 
              for (j = 0; j < path_length; j++) {
                path_dir[j] = dirs[((k + j) + 1) % length];
              }
            if (dirs[k] > 3) 
              for (j = 0; j < path_length; j++) {
                path_dir[(path_length - 1) - j] = (7 - dirs[((k + j) + 1) % length]);
              }
/**if(dir==XUP)printf("X_UPDATE PATH: "); printpath( path_dir, path_length );**/
            path_product(path_dir,path_length);
/* We took the path in the other direction from our
			old convention in order to get it to end up
			"at our site", so now take adjoint */
/* then compute "single_action" contribution to
			staple */
            for (((i = 0) , (st = lattice)); i < sites_on_node; (i++ , st++)) {
              su3_adjoint(&st -> tempmat1,&tmat1);
/* first we compute the fundamental term */
              new_term = loop_coeff[iloop][0];
/* now we add in the higher representations */
              if (0) {
                if (this_node == 0) 
                  printf("WARNING: THIS CODE IS NOT TESTED\n");
                exit(0);
                act2 = 1.0;
                action = (3.0 - (realtrace_su3(((st -> link) + dir),&tmat1)));
                for (j = 1; j < 1; j++) {
                  act2 *= action;
                  new_term += ((loop_coeff[iloop][j] * act2) * ((float )(j + 1)));
                }
/* end if NREPS > 1 */
              }
              scalar_mult_add_su3_matrix(&st -> staple,&tmat1,new_term,&st -> staple);
            }
            ncount++;
/* k (location in path) */
          }
/* ln */
      }
/* iloop */
    }
/* Now multiply the staple sum by the link, then update momentum */
    for (((i = 0) , (st = lattice)); i < sites_on_node; (i++ , st++)) {
      mult_su3_na(((st -> link) + dir),&st -> staple,&tmat1);
      momentum = ((anti_hermitmat *)(((char *)st) + mom_off));
      uncompress_anti_hermitian((momentum + dir),&tmat2);
      scalar_mult_sub_su3_matrix(&tmat2,&tmat1,eb3,&st -> staple);
      make_anti_hermitian(&st -> staple,(momentum + dir));
    }
/* dir loop */
  }
#ifdef GFTIME
  dtime += dclock();
  if (this_node == 0) 
    printf("GFTIME:   time = %e (Symanzik1) mflops = %e\n",dtime,((nflop * volume) / ((1e6 * dtime) * (numnodes()))));
#endif
/* imp_gauge_force.c */
}
/* Measure gauge observables:
    Loops in action (time and space directions treated differently)
    Polyakov loop
*/

void g_measure()
{
  double ss_plaquette;
  double st_plaquette;
  complex p_loop;
  register int i;
  register site *s;
  complex trace;
  double average[1UL];
  double action;
  double act2;
  double total_action;
  int length;
/* these are for loop_table  */
  int ln;
  int iloop;
  int rep;
/* KS and BC minus signs should be out for this routine */
  d_plaquette(&ss_plaquette,&st_plaquette);
  if (this_node == 0) 
    printf("PLAQ:\t%f\t%f\n",ss_plaquette,st_plaquette);
  p_loop = ploop();
  if (this_node == 0) 
    printf("P_LOOP:\t%e\t%e\n",p_loop.real,p_loop.imag);
/* gauge action, all loops that contribute */
  total_action = 0.0;
  for (iloop = 0; iloop < 3; iloop++) {
    length = loop_length[iloop];
/* loop over rotations and reflections */
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      path_product(loop_table[iloop][ln],length);
      for (rep = 0; rep < 1; rep++) 
        average[rep] = 0.0;
      for (((i = 0) , (s = lattice)); i < sites_on_node; (i++ , s++)) {
        trace = trace_su3(&s -> tempmat1);
        average[0] += ((double )trace.real);
        action = (3.0 - ((double )trace.real));
        total_action += (((double )loop_coeff[iloop][0]) * action);
/* need the "3 -" for higher characters */
        act2 = action;
        for (rep = 1; rep < 1; rep++) {
          act2 *= action;
          average[rep] += act2;
          total_action += (((double )loop_coeff[iloop][rep]) * act2);
/* reps */
        }
/* sites */
      }
      g_vecdoublesum(average,1);
/* dump the loop */
/*
	    node0_printf("G_LOOP:  %d  %d  %d   ",iloop,ln,length);
	    for(rep=0;rep<NREPS;rep++)node0_printf("\t%e",average[rep]/volume);
	    node0_printf("\t( ");
	    for(i=0;i<length;i++)node0_printf("%d ",loop_table[iloop][ln][i]);
	    node0_printf(" )\n");
	    */
/* ln */
    }
/* iloop */
  }
  g_doublesum(&total_action);
/* node0_printf("GACTION: %e\n",total_action/volume); */
/**node0_printf("CHECK:   %e   %e\n",total_action,imp_gauge_action() );**/
  if (this_node == 0) 
    fflush(stdout);
}

void printpath(int *path,int length)
{
  register int i;
  if (this_node == 0) 
    printf("\t( ");
  for (i = 0; i < length; i++) 
    if (this_node == 0) 
      printf("%d ",path[i]);
  if (this_node == 0) 
    printf(",  L = %d )\n",length);
}
#ifdef N_SUBL32
/*** code from symanzik_sl32/dsdu_qhb.c  -- compute the staple ***/
/* This is a version for extended actions where 32 sublattices are
   needed to make the links independent. */
/* U.M. Heller August 1997 */
/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */
#include <assert.h>
/* This procedure designed only for NREPS = 1 */
/* set up dirs.  we are looking at loop starting in "XUP"
	       direction, rotate so it starts in "dir" direction. */
/* generalized "staple" */
/* The path starts at the forward end of the link */
/* check for links in direction of link to be updated.
	       Note the direction of the path - opposite the link. */
/* We took the path in the other direction from our old
		   convention in order to get it to end up "at our site".
		   So now take adjoint */
/* k (location in path) */
/* ln */
/* iloop */
/* dsdu_qhb */
#endif /* N_SUBL32 */
