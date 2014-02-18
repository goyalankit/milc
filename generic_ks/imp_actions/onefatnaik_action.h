#ifndef _ONEFATNAIK_ACTION_H
#define _ONEFATNAIK_ACTION_H

#include "../include/dirs.h"

    /* The one-link plus Naik plus staple(fat) action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define MAX_BASIC_PATHS 3
#define MAX_LENGTH 3
#define MAX_NUM 64


    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR },
    { XUP, XUP, XUP },
    { YUP, XUP, YDOWN }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static float path_coeff[MAX_BASIC_PATHS] = {
        (9.0/8.0)/4.0,	/* one link */
       -(9.0/8.0)/27.0,     /* three link (Naik) */
       -(9.0/8.0)/8.0      /* staple(fat) */
    };
    static char quark_action_description[] =
      "\"one-link + Naik + staple(fat) action\"";
#endif /* _ONEFATNAIK_ACTION_H */
