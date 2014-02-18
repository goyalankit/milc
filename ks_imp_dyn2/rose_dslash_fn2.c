/******* dslash_fn2.c - dslash for improved KS fermions T3E version ****/
/* MIMD version 6 */
/* Kogut-Susskind fermions -- improved.  This version for "fat plus
   Naik" quark action.  Connection to nearest neighbors stored in
   fatlink and to third nearest neighbors in longlink */
/* With DSLASH_TMP_LINKS, assumes that the gauge links have been
   prestored in t_fatlinks and t_longlinks.  Otherwise, takes the
   fatlinks and longlinks from the site structure. */
/* This version waits for gathers from both positive and negative
   directions before computing, thereby combining two lattice loops in
   an attempt to gain prefetching time for sub_four_su3_vecs */
/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* C. DeTar 9/29/01 Standardized prefetching and synced the versions */
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define INDEX_3RD(dir) (dir - 8)      /* this gives the 'normal' direction */
/* Temporary work space for dslash_fn_on_temp_special */
static su3_vector *temp[9UL];
/* Flag indicating if temp is allocated               */
static int temp_not_allocated = 1;

void cleanup_gathers(msg_tag *tags1[],msg_tag *tags2[])
{
  int i;
  for (i = 0; i <= 3; i++) {
    cleanup_gather(tags1[i]);
    cleanup_gather(tags1[7 - i]);
    cleanup_gather(tags2[i]);
    cleanup_gather(tags2[7 - i]);
  }
  for (i = 8; i <= 11; i++) {
    cleanup_gather(tags1[i]);
    cleanup_gather(tags1[23 - i]);
    cleanup_gather(tags2[i]);
    cleanup_gather(tags2[23 - i]);
  }
}

void cleanup_dslash_temps()
{
  register int i;
  if (!(temp_not_allocated != 0)) 
    for (i = 0; i < 9; i++) {
      free(temp[i]);
    }
  temp_not_allocated = 1;
}
/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.  Use "fatlinks" for one link transport,
   "longlinks" for three link transport. */

void dslash_fn(field_offset src,field_offset dest,int parity)
{
  register int i;
  register site *s;
  register int dir;
  register int otherparity;
  register su3_matrix *fat4;
  register su3_matrix *long4;
  msg_tag *tag[16UL];
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
  switch(parity){
    case 0x02:
{
      otherparity = 1;
      break; 
    }
    case 0x01:
{
      otherparity = 2;
      break; 
    }
    case 0x03:
{
      otherparity = 3;
      break; 
    }
  }
/* Start gathers from positive directions */
/* And start the 3-step gather too */
  for (dir = 0; dir <= 3; dir++) {
    tag[dir] = start_gather(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir]);
    tag[dir + 8] = start_gather(src,(sizeof(su3_vector )),(dir + 8),parity,gen_pt[dir + 8]);
  }
/* Multiply by adjoint matrix at other sites */
/* Use fat link for single link transport */
{
    register int loopend;
    loopend = ((otherparity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_adj_su3_mat_vec_4dir(fat4,((su3_vector *)(((char *)s) + src)),(s -> tempvec));
/* multiply by 3-link matrices too */
      mult_adj_su3_mat_vec_4dir(long4,((su3_vector *)(((char *)s) + src)),(s -> templongvec));
    }
  }
/* Start gathers from negative directions */
  for (dir = 0; dir <= 3; dir++) {
    tag[7 - dir] = start_gather(((field_offset )(((char *)(lattice[0].tempvec + dir)) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir]);
  }
/* Start 3-neighbour gathers from negative directions */
  for (dir = 8; dir <= 11; dir++) {
    tag[23 - dir] = start_gather(((field_offset )(((char *)(lattice[0].templongvec + (dir - 8))) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir]);
  }
/* Wait gathers from positive directions, multiply by matrix and
	accumulate */
/* wait for the 3-neighbours from positive directions, multiply */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[dir]);
    wait_gather(tag[dir + 8]);
  }
/* Wait gathers from negative directions, accumulate (negative) */
/* and the same for the negative 3-rd neighbours */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    wait_gather(tag[23 - dir]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_su3_mat_vec_sum_4dir(fat4,((su3_vector *)gen_pt[0][i]),((su3_vector *)gen_pt[1][i]),((su3_vector *)gen_pt[2][i]),((su3_vector *)gen_pt[3][i]),((su3_vector *)(((char *)s) + dest)));
      mult_su3_mat_vec_sum_4dir(long4,((su3_vector *)gen_pt[8][i]),((su3_vector *)gen_pt[9][i]),((su3_vector *)gen_pt[10][i]),((su3_vector *)gen_pt[11][i]),&s -> templongv1);
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;;
      }
      sub_four_su3_vecs(((su3_vector *)(((char *)s) + dest)),((su3_vector *)gen_pt[7][i]),((su3_vector *)gen_pt[6][i]),((su3_vector *)gen_pt[5][i]),((su3_vector *)gen_pt[4][i]));
      sub_four_su3_vecs(&s -> templongv1,((su3_vector *)gen_pt[15][i]),((su3_vector *)gen_pt[14][i]),((su3_vector *)gen_pt[13][i]),((su3_vector *)gen_pt[12][i]));
/* Now need to add these things together */
      add_su3_vector(((su3_vector *)(((char *)s) + dest)),&s -> templongv1,((su3_vector *)(((char *)s) + dest)));
    }
  }
/* free up the buffers */
  for (dir = 0; dir <= 3; dir++) {
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[23 - dir]);
  }
}
/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather, otherwise use restart_gather. 
  The calling program must clean up the gathers! */

void dslash_fn_special(field_offset src,field_offset dest,int parity,msg_tag **tag,int start)
{
  register int i;
  register site *s;
  register int dir;
  register int otherparity;
  register su3_matrix *fat4;
  register su3_matrix *long4;
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
  switch(parity){
    case 0x02:
{
      otherparity = 1;
      break; 
    }
    case 0x01:
{
      otherparity = 2;
      break; 
    }
    case 0x03:
{
      otherparity = 3;
      break; 
    }
  }
/* Start gathers from positive directions */
  for (dir = 0; dir <= 3; dir++) {
/**printf("dslash_special: up gathers, start=%d\n",start);**/
    if (start == 1) 
      tag[dir] = start_gather(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir]);
    else 
      restart_gather(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir],tag[dir]);
  }
/* and start the 3rd neighbor gather */
  for (dir = 8; dir <= 11; dir++) {
    if (start == 1) 
      tag[dir] = start_gather(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir]);
    else 
      restart_gather(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir],tag[dir]);
  }
/* Multiply by adjoint matrix at other sites */
{
    register int loopend;
    loopend = ((otherparity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_adj_su3_mat_vec_4dir(fat4,((su3_vector *)(((char *)s) + src)),(s -> tempvec));
/* multiply by 3-link matrices too */
      mult_adj_su3_mat_vec_4dir(long4,((su3_vector *)(((char *)s) + src)),(s -> templongvec));
    }
  }
/* Start gathers from negative directions */
  for (dir = 0; dir <= 3; dir++) {
/**printf("dslash_special: down gathers, start=%d\n",start);**/
    if (start == 1) 
      tag[7 - dir] = start_gather(((field_offset )(((char *)(lattice[0].tempvec + dir)) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir]);
    else 
      restart_gather(((field_offset )(((char *)(lattice[0].tempvec + dir)) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir],tag[7 - dir]);
  }
/* and 3rd neighbours */
  for (dir = 8; dir <= 11; dir++) {
/**printf("dslash_special: down gathers, start=%d\n",start);**/
    if (start == 1) 
      tag[23 - dir] = start_gather(((field_offset )(((char *)(lattice[0].templongvec + (dir - 8))) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir]);
    else 
      restart_gather(((field_offset )(((char *)(lattice[0].templongvec + (dir - 8))) - ((char *)(lattice + 0)))),(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir],tag[23 - dir]);
  }
/* Wait gathers from positive directions, multiply by matrix and
	accumulate */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[dir]);
  }
/* wait for the 3-neighbours from positive directions, multiply */
  for (dir = 8; dir <= 11; dir++) {
    wait_gather(tag[dir]);
  }
/* Wait gathers from negative directions, accumulate (negative) */
/* and the same for the negative 3-rd neighbours */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    wait_gather(tag[23 - dir]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_su3_mat_vec_sum_4dir(fat4,((su3_vector *)gen_pt[0][i]),((su3_vector *)gen_pt[1][i]),((su3_vector *)gen_pt[2][i]),((su3_vector *)gen_pt[3][i]),((su3_vector *)(((char *)s) + dest)));
      mult_su3_mat_vec_sum_4dir(long4,((su3_vector *)gen_pt[8][i]),((su3_vector *)gen_pt[9][i]),((su3_vector *)gen_pt[10][i]),((su3_vector *)gen_pt[11][i]),&s -> templongv1);
      sub_four_su3_vecs(((su3_vector *)(((char *)s) + dest)),((su3_vector *)gen_pt[7][i]),((su3_vector *)gen_pt[6][i]),((su3_vector *)gen_pt[5][i]),((su3_vector *)gen_pt[4][i]));
      sub_four_su3_vecs(&s -> templongv1,((su3_vector *)gen_pt[15][i]),((su3_vector *)gen_pt[14][i]),((su3_vector *)gen_pt[13][i]),((su3_vector *)gen_pt[12][i]));
/*** Now need to add these things together ***/
      add_su3_vector(((su3_vector *)(((char *)s) + dest)),&s -> templongv1,((su3_vector *)(((char *)s) + dest)));
    }
  }
}

void dslash_fn_on_temp(su3_vector *src,su3_vector *dest,int parity)
{
  register int i;
  register site *s;
  register int dir;
  register int otherparity;
  msg_tag *tag[16UL];
  su3_vector *tempvec[4UL];
  su3_vector *templongvec[4UL];
  su3_vector *templongv1;
  register su3_matrix *fat4;
  register su3_matrix *long4;
  for (dir = 0; dir <= 3; dir++) {
    tempvec[dir] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    templongvec[dir] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
  }
  templongv1 = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
  switch(parity){
    case 0x02:
{
      otherparity = 1;
      break; 
    }
    case 0x01:
{
      otherparity = 2;
      break; 
    }
    case 0x03:
{
      otherparity = 3;
      break; 
    }
  }
/* Start gathers from positive directions */
/* And start the 3-step gather too */
  for (dir = 0; dir <= 3; dir++) {
    tag[dir] = start_gather_from_temp(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir]);
    tag[dir + 8] = start_gather_from_temp(src,(sizeof(su3_vector )),(dir + 8),parity,gen_pt[dir + 8]);
  }
/* Multiply by adjoint matrix at other sites */
/* Use fat link for single link transport */
{
    register int loopend;
    loopend = ((otherparity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_adj_su3_mat_4vec(fat4,(src + i),(tempvec[0] + i),(tempvec[1] + i),(tempvec[2] + i),(tempvec[3] + i));
/* multiply by 3-link matrices too */
      mult_adj_su3_mat_4vec(long4,(src + i),(templongvec[0] + i),(templongvec[1] + i),(templongvec[2] + i),(templongvec[3] + i));
    }
  }
/* Start gathers from negative directions */
  for (dir = 0; dir <= 3; dir++) {
    tag[7 - dir] = start_gather_from_temp(tempvec[dir],(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir]);
  }
/* Start 3-neighbour gathers from negative directions */
  for (dir = 8; dir <= 11; dir++) {
    tag[23 - dir] = start_gather_from_temp(templongvec[dir - 8],(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir]);
  }
/* Wait gathers from positive directions, multiply by matrix and
	accumulate */
/* wait for the 3-neighbours from positive directions, multiply */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[dir]);
    wait_gather(tag[dir + 8]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_su3_mat_vec_sum_4dir(fat4,((su3_vector *)gen_pt[0][i]),((su3_vector *)gen_pt[1][i]),((su3_vector *)gen_pt[2][i]),((su3_vector *)gen_pt[3][i]),(dest + i));
      mult_su3_mat_vec_sum_4dir(long4,((su3_vector *)gen_pt[8][i]),((su3_vector *)gen_pt[9][i]),((su3_vector *)gen_pt[10][i]),((su3_vector *)gen_pt[11][i]),(templongv1 + i));
    }
  }
/* Wait gathers from negative directions, accumulate (negative) */
/* and the same for the negative 3-rd neighbours */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    wait_gather(tag[23 - dir]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {;;;
      }
      sub_four_su3_vecs((dest + i),((su3_vector *)gen_pt[7][i]),((su3_vector *)gen_pt[6][i]),((su3_vector *)gen_pt[5][i]),((su3_vector *)gen_pt[4][i]));
      sub_four_su3_vecs((templongv1 + i),((su3_vector *)gen_pt[15][i]),((su3_vector *)gen_pt[14][i]),((su3_vector *)gen_pt[13][i]),((su3_vector *)gen_pt[12][i]));
/* Now need to add these things together */
      add_su3_vector((dest + i),(templongv1 + i),(dest + i));
    }
  }
/* free up the buffers */
  for (dir = 0; dir <= 3; dir++) {
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[23 - dir]);
  }
  for (dir = 0; dir <= 3; dir++) {
    free(tempvec[dir]);
    free(templongvec[dir]);
  }
  free(templongv1);
}
/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Next to last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather, otherwise use restart_gather. 
  The calling program must clean up the gathers and temps! */

void dslash_fn_on_temp_special(su3_vector *src,su3_vector *dest,int parity,msg_tag **tag,int start)
{
  register int i;
  register site *s;
  register int dir;
  register int otherparity;
  register su3_matrix *fat4;
  register su3_matrix *long4;
/* allocate temporary work space only if not already allocated */
  if (temp_not_allocated != 0) {
    for (dir = 0; dir <= 3; dir++) {
      temp[dir] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
      temp[dir + 4] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    }
    temp[8] = ((su3_vector *)(malloc((sites_on_node * sizeof(su3_vector )))));
    temp_not_allocated = 0;
  }
/* load fatlinks and longlinks */
  if (!(valid_longlinks != 0)) 
    load_longlinks();
  if (!(valid_fatlinks != 0)) 
    load_fatlinks();
  switch(parity){
    case 0x02:
{
      otherparity = 1;
      break; 
    }
    case 0x01:
{
      otherparity = 2;
      break; 
    }
    case 0x03:
{
      otherparity = 3;
      break; 
    }
  }
/* Start gathers from positive directions */
/* And start the 3-step gather too */
  for (dir = 0; dir <= 3; dir++) {
    if (start == 1) {
      tag[dir] = start_gather_from_temp(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir]);
      tag[dir + 8] = start_gather_from_temp(src,(sizeof(su3_vector )),(dir + 8),parity,gen_pt[dir + 8]);
    }
    else {
      restart_gather_from_temp(src,(sizeof(su3_vector )),dir,parity,gen_pt[dir],tag[dir]);
      restart_gather_from_temp(src,(sizeof(su3_vector )),(dir + 8),parity,gen_pt[dir + 8],tag[dir + 8]);
    }
  }
/* Multiply by adjoint matrix at other sites */
/* Use fat link for single link transport */
{
    register int loopend;
    loopend = ((otherparity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((otherparity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_adj_su3_mat_4vec(fat4,(src + i),(temp[0] + i),(temp[1] + i),(temp[2] + i),(temp[3] + i));
/* multiply by 3-link matrices too */
      mult_adj_su3_mat_4vec(long4,(src + i),(temp[4] + i),(temp[5] + i),(temp[6] + i),(temp[7] + i));
    }
  }
/* Start gathers from negative directions */
  for (dir = 0; dir <= 3; dir++) {
    if (start == 1) 
      tag[7 - dir] = start_gather_from_temp(temp[dir],(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir]);
    else 
      restart_gather_from_temp(temp[dir],(sizeof(su3_vector )),(7 - dir),parity,gen_pt[7 - dir],tag[7 - dir]);
  }
/* Start 3-neighbour gathers from negative directions */
  for (dir = 8; dir <= 11; dir++) {
    if (start == 1) 
      tag[23 - dir] = start_gather_from_temp(temp[(dir - 8) + 4],(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir]);
    else 
      restart_gather_from_temp(temp[(dir - 8) + 4],(sizeof(su3_vector )),(23 - dir),parity,gen_pt[23 - dir],tag[23 - dir]);
  }
/* Wait gathers from positive directions, multiply by matrix and
	accumulate */
/* wait for the 3-neighbours from positive directions, multiply */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[dir]);
    wait_gather(tag[dir + 8]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {
#ifdef DSLASH_TMP_LINKS
        fat4 = (t_fatlink + (4 * (i + 1)));
        long4 = (t_longlink + (4 * (i + 1)));
#else
#endif
;;;;
      }
#ifdef DSLASH_TMP_LINKS
      fat4 = (t_fatlink + (4 * i));
      long4 = (t_longlink + (4 * i));
#else
#endif
      mult_su3_mat_vec_sum_4dir(fat4,((su3_vector *)gen_pt[0][i]),((su3_vector *)gen_pt[1][i]),((su3_vector *)gen_pt[2][i]),((su3_vector *)gen_pt[3][i]),(dest + i));
      mult_su3_mat_vec_sum_4dir(long4,((su3_vector *)gen_pt[8][i]),((su3_vector *)gen_pt[9][i]),((su3_vector *)gen_pt[10][i]),((su3_vector *)gen_pt[11][i]),(temp[8] + i));
    }
  }
/* Wait gathers from negative directions, accumulate (negative) */
/* and the same for the negative 3-rd neighbours */
  for (dir = 0; dir <= 3; dir++) {
    wait_gather(tag[7 - dir]);
  }
  for (dir = 8; dir <= 11; dir++) {
    wait_gather(tag[23 - dir]);
  }
{
    register int loopend;
    loopend = ((parity == 2)?even_sites_on_node : sites_on_node);
    for (((i = ((parity == 1)?even_sites_on_node : 0)) , (s = (lattice + i))); i < loopend; (i++ , s++)) {
      if (i < (loopend - 1)) {;;
      }
      sub_four_su3_vecs((dest + i),((su3_vector *)gen_pt[7][i]),((su3_vector *)gen_pt[6][i]),((su3_vector *)gen_pt[5][i]),((su3_vector *)gen_pt[4][i]));
      sub_four_su3_vecs((temp[8] + i),((su3_vector *)gen_pt[15][i]),((su3_vector *)gen_pt[14][i]),((su3_vector *)gen_pt[13][i]),((su3_vector *)gen_pt[12][i]));
/* Now need to add these things together */
      add_su3_vector((dest + i),(temp[8] + i),(dest + i));
    }
  }
}
