/********************** io_helpers.c **********************************/
/* MIMD version 6 */
/* DT 8/97 
     General purpose high level routines, to be used by any application
     that wants them.
*/
#include "generic_includes.h"
#include "../include/io_lat.h"
extern FILE *fpInput;
/* save a lattice in any of the formats:
    SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_CHECKPOINT
*/

gauge_file *save_lattice(int flag,char *filename)
{
  double dtime;
  gauge_file *gf;
  double ssplaq;
  double stplaq;
  dtime = -dclock();
  switch(flag){
    case 20:
{
      gf = ((gauge_file *)((void *)0));
      break; 
    }
    case 21:
{
      gf = save_ascii(filename);
      break; 
    }
    case 22:
{
      gf = save_serial(filename);
      break; 
    }
    case 28:
{
      gf = save_parallel(filename);
      break; 
    }
    case 23:
{
      gf = save_checkpoint(filename);
      break; 
    }
    case 30:
{
      gf = save_serial_archive(filename);
      break; 
    }
    case 31:
{
      gf = save_parallel_archive(filename);
      break; 
    }
    default:
{
      printf("\nsave_lattice: ERROR: unknown type for saving lattice\n");
      terminate(1);
    }
  }
  dtime += dclock();
  if (flag != 20) 
    if (this_node == 0) 
      printf("Time to save = %e\n",dtime);
  d_plaquette(&ssplaq,&stplaq);
  if (this_node == 0) 
    printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);
  return gf;
}
/* reload a lattice in any of the formats, or cold lattice, or keep
     current lattice:
    FRESH, CONTINUE,
    RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL
*/

gauge_file *reload_lattice(int flag,char *filename)
{
  double dtime;
  gauge_file *gf;
  double ssplaq;
  double stplaq;
  float max_deviation;
  void coldlat();
  dtime = -dclock();
  switch(flag){
/* do nothing */
    case 10:
{
      gf = ((gauge_file *)((void *)0));
      break; 
    }
/* cold lattice */
    case 11:
{
      coldlat();
      gf = ((gauge_file *)((void *)0));
      break; 
    }
/* read Ascii lattice */
    case 12:
{
      gf = restore_ascii(filename);
      break; 
    }
/* read binary lattice serially */
    case 13:
{
      gf = restore_serial(filename);
      break; 
    }
/* read binary lattice in parallel */
    case 19:
{
      gf = restore_parallel(filename);
      break; 
    }
    default:
{
      if (this_node == 0) 
        printf("reload_lattice: Bad startflag %d\n",flag);
      terminate(1);
    }
  }
  dtime += dclock();
  if ((flag != 11) && (flag != 10)) 
    if (this_node == 0) 
      printf("Time to reload gauge configuration = %e\n",dtime);
#ifdef SCHROED_FUN
#endif
  d_plaquette(&ssplaq,&stplaq);
  if (this_node == 0) {
    printf("CHECK PLAQ: %e %e\n",ssplaq,stplaq);
    fflush(stdout);
  }
  dtime = -dclock();
  max_deviation = check_unitarity();
  g_floatmax(&max_deviation);
  dtime += dclock();
  if (this_node == 0) 
    printf("Unitarity checked.  Max deviation %.2e\n",max_deviation);
  fflush(stdout);
  if (this_node == 0) 
    printf("Time to check unitarity = %e\n",dtime);
  return gf;
}
/* find out what kind of starting lattice to use, and lattice name if
   necessary.  This routine is only called by node 0.
*/

int ask_starting_lattice(int prompt,int *flag,char *filename)
{
  char savebuf[256UL];
  int status;
  if (prompt != 0) 
    printf("enter \'continue\', \'fresh\', \'reload_ascii\', \'reload_serial\', or \'reload_parallel\'\n");
  status = fscanf(fpInput,"%s",savebuf);
  if (status != 1) {
    printf("\nask_starting_lattice: ERROR IN INPUT: starting lattice command \"%s\" is invalid\n",savebuf);
    return 1;
  }
  printf("%s ",savebuf);
  if (strcmp("fresh",savebuf) == 0) {
     *flag = 11;
    printf("\n");
  }
  else if (strcmp("continue",savebuf) == 0) {
     *flag = 10;
    printf("\n");
  }
  else if (strcmp("reload_ascii",savebuf) == 0) {
     *flag = 12;
  }
  else if (strcmp("reload_serial",savebuf) == 0) {
     *flag = 13;
  }
  else if (strcmp("reload_parallel",savebuf) == 0) {
     *flag = 19;
  }
  else {
    printf("\nask_starting_lattice: ERROR IN INPUT: lattice_command \"%s\" is invalid\n",savebuf);
    return 1;
  }
/*read name of file and load it */
  if (( *flag != 11) && ( *flag != 10)) {
    if (prompt != 0) 
      printf("enter name of file containing lattice\n");
    status = fscanf(fpInput,"%s",filename);
    if (status != 1) {
      printf("\nask_starting_lattice: ERROR IN INPUT: error reading file name\n");
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
}
/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/

int ask_ending_lattice(int prompt,int *flag,char *filename)
{
  char savebuf[256UL];
  int status;
  if (prompt != 0) 
    printf("\'forget\' lattice at end,  \'save_ascii\', \'save_serial\', \'save_parallel\', \'save_checkpoint\', \'save_serial_archive\', or \'save_parallel_archive\'\n");
  status = fscanf(fpInput,"%s",savebuf);
  if (status != 1) {
    printf("\nask_ending_lattice: ERROR IN INPUT: error reading ending lattice command\n");
    return 1;
  }
  printf("%s ",savebuf);
  if (strcmp("save_ascii",savebuf) == 0) {
     *flag = 21;
  }
  else if (strcmp("save_serial",savebuf) == 0) {
     *flag = 22;
  }
  else if (strcmp("save_parallel",savebuf) == 0) {
     *flag = 28;
  }
  else if (strcmp("save_checkpoint",savebuf) == 0) {
     *flag = 23;
  }
  else if (strcmp("save_serial_archive",savebuf) == 0) {
     *flag = 30;
  }
  else if (strcmp("save_parallel_archive",savebuf) == 0) {
     *flag = 31;
  }
  else if (strcmp("forget",savebuf) == 0) {
     *flag = 20;
    printf("\n");
  }
  else {
    printf("\nask_ending_lattice: ERROR IN INPUT: \"%s\" is not a save lattice command\n",savebuf);
    return 1;
  }
  if ( *flag != 20) {
    if (prompt != 0) 
      printf("enter filename\n");
    status = fscanf(fpInput,"%s",filename);
    if (status != 1) {
      printf("\nask_ending_lattice: ERROR IN INPUT: error reading filename\n");
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
}

void coldlat()
{
/* sets link matrices to unit matrices */
  register int i;
  register int j;
  register int k;
  register int dir;
  register site *sit;
  for (((i = 0) , (sit = lattice)); i < sites_on_node; (i++ , sit++)) {
    for (dir = 0; dir <= 3; dir++) {
      for (j = 0; j < 3; j++) {
        for (k = 0; k < 3; k++) {
          if (j != k) {
            (sit -> link)[dir].e[j][k] = cmplx(0.0,0.0);
          }
          else {
            (sit -> link)[dir].e[j][k] = cmplx(1.0,0.0);
          }
        }
      }
    }
  }
  if (this_node == 0) 
    printf("unit gauge configuration loaded\n");
}

void funnylat()
{
/* sets link matrices to funny matrices for debugging */
  register int i;
  register int j;
  register int k;
  register int dir;
  register site *sit;
  for (((i = 0) , (sit = lattice)); i < sites_on_node; (i++ , sit++)) {
    for (dir = 0; dir <= 3; dir++) {
      for (j = 0; j < 3; ++j) {
        for (k = 0; k < 3; ++k) {
          (sit -> link)[dir].e[j][k] = cmplx(0.0,0.0);
        }
      }
      (sit -> link)[dir].e[0][0].real = dir;
      (sit -> link)[dir].e[1][1].real = (10 * (sit -> x));
      (sit -> link)[dir].e[2][2].real = (100 * (sit -> y));
      (sit -> link)[dir].e[0][0].imag = dir;
      (sit -> link)[dir].e[1][1].imag = (10 * (sit -> z));
      (sit -> link)[dir].e[2][2].imag = (100 * (sit -> t));
    }
  }
}
/* get_f is used to get a floating point number.  If prompt is non-zero,
it will prompt for the input value with the variable_name_string.  If
prompt is zero, it will require that variable_name_string precede the
input value.  get_i gets an integer.
get_i and get_f return the values, and exit on error */

int get_f(int prompt,char *variable_name_string,float *value)
{
  int s;
  char checkname[80UL];
  if (prompt != 0) {
    s = 0;
    while(s != 1){
      printf("enter %s ",variable_name_string);
      fscanf(fpInput,"%s",checkname);
      s = sscanf(checkname,"%e",value);
      if (s == 1) 
        printf("%s %g\n",variable_name_string,( *value));
      else 
        printf("Data format error.\n");
    }
  }
  else {
    s = fscanf(fpInput,"%s",checkname);
    if (s == -1) {
      printf("get_f: EOF on STDIN while expecting %s.\n",variable_name_string);
      return 1;
    }
    else if (s == 0) {
      printf("get_f: Format error looking for %s\n",variable_name_string);
      return 1;
    }
    else if (strcmp(checkname,variable_name_string) != 0) {
      printf("\nget_f: ERROR IN INPUT: expected %s but found %s\n",variable_name_string,checkname);
      return 1;
    }
    printf("%s ",variable_name_string);
    s = fscanf(fpInput,"%e",value);
    if (s == -1) {
      printf("\nget_f: Expecting value for %s but found EOF.\n",variable_name_string);
      return 1;
    }
    else if (s == 0) {
      printf("\nget_f: Format error reading value for %s\n",variable_name_string);
      return 1;
    }
    printf("%g\n",( *value));
  }
  return 0;
}

int get_i(int prompt,char *variable_name_string,int *value)
{
  int s;
  char checkname[80UL];
  if (prompt != 0) {
    s = 0;
    while(s != 1){
      printf("enter %s ",variable_name_string);
      fscanf(fpInput,"%s",checkname);
      s = sscanf(checkname,"%d",value);
      if (s == 1) 
        printf("%s %d\n",variable_name_string, *value);
      else 
        printf("Data format error.\n");
    }
  }
  else {
    s = fscanf(fpInput,"%s",checkname);
    if (s == -1) {
/*
	  printf("get_i: EOF on STDIN while expecting %s.\n",
	       variable_name_string);
	*/
      return 1;
    }
    else if (s == 0) {
      printf("\nget_i: Format error looking for %s\n",variable_name_string);
      return 1;
    }
    else if (strcmp(checkname,variable_name_string) != 0) {
      printf("\nget_i: ERROR IN INPUT: expected %s but found %s\n",variable_name_string,checkname);
      return 1;
    }
    printf("%s ",variable_name_string);
    s = fscanf(fpInput,"%d",value);
    if (s == -1) {
      printf("\nget_i: Expecting value for %s but found EOF.\n",variable_name_string);
      return 1;
    }
    else if (s == 0) {
      printf("\nget_i: Format error reading value for %s\n",variable_name_string);
      return 1;
    }
    printf("%d\n", *value);
  }
  return 0;
}
/* Read a single word as a string */

int get_s(int prompt,char *variable_name_string,char *value)
{
  int s;
  char checkname[80UL];
  if (prompt != 0) {
    s = 0;
    while(s != 1){
      printf("enter %s ",variable_name_string);
      s = fscanf(fpInput,"%s",value);
      if (s == 1) 
        printf("%s %s\n",variable_name_string,value);
      else 
        printf("Data format error.\n");
    }
  }
  else {
    s = fscanf(fpInput,"%s",checkname);
    if (s == -1) {
      printf("get_s: EOF on STDIN while expecting %s.\n",variable_name_string);
      return 1;
    }
    else if (s == 0) {
      printf("\nget_s: Format error looking for %s\n",variable_name_string);
      return 1;
    }
    else if (strcmp(checkname,variable_name_string) != 0) {
      printf("\nget_s: ERROR IN INPUT: expected %s but found %s\n",variable_name_string,checkname);
      return 1;
    }
    printf("%s ",variable_name_string);
    s = fscanf(fpInput,"%s",value);
    if (s == -1) {
      printf("\nget_s: Expecting value for %s but found EOF.\n",variable_name_string);
      return 1;
    }
    else if (s == 0) {
      printf("\nget_s: Format error reading value for %s\n",variable_name_string);
      return 1;
    }
    printf("%s\n",value);
  }
  return 0;
}
/* get_prompt gets the initial value of prompt */
/* 0 for reading from file, 1 prompts for input from terminal */
/* should be called only by node 0 */
/* return 0 if sucessful, 1 if failure */

int get_prompt(int *prompt)
{
  char initial_prompt[80UL];
   *prompt = -1;
  printf("type 0 for no prompts  or 1 for prompts\n");
  fscanf(fpInput,"%s",initial_prompt);
  if (strcmp(initial_prompt,"prompt") == 0) {
    fscanf(fpInput,"%d",prompt);
  }
  else if (strcmp(initial_prompt,"0") == 0) 
     *prompt = 0;
  else if (strcmp(initial_prompt,"1") == 0) 
     *prompt = 1;
  if (( *prompt == 0) || ( *prompt == 1)) 
    return 0;
  else {
    printf("\nget_prompt: ERROR IN INPUT: initial prompt\n");
    return 1;
  }
}
