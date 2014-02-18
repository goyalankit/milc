#ifndef _GENERIC_PG_H
#define _GENERIC_PG_H
/************************ generic_pg.h **********************************
*									*
*  Macros and declarations for generic_pg routines                      *
*  This header is for codes that call generic_pg routines               *
*  MIMD version 6 							*
*									*
*/

#include "../include/macros.h"

int update();
void update_h(float eps);
void update_u(float eps);
void relax(int NumStp);
void monte(int NumStp);
void dsdu_qhb(int dir1, int parity);
double d_action();
void gauge_field_copy(field_offset src, field_offset dest);
#endif	/* _GENERIC_PG_H */
