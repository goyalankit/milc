#ifndef _FAT7TAD_ACTION_H
#define _FAT7TAD_ACTION_H

#include "../include/dirs.h"

    /* The fat link action with five link paths designed to minimize
	couplings at momentum pi in any direction.  */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define MAX_BASIC_PATHS 4
#define MAX_LENGTH 7
#define MAX_NUM 632
#define TADPOLE_IMPROVE	/* use tadpole improvement in quark action */

    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,5,7};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static float path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0/8.0),                  /* one link */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
    };
    static char quark_action_description[] =
	"\"Fat-7: seven link paths, couplings(pi)=0, tadpole weights\"";

#endif /* _FAT7TAD_ACTION_H */
