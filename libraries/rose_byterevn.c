/******************************** byterevn.c ***************************/
/* MIMD version 6 */
/* WARNING - MUST BE COMPILED WITH APPROPRIATE SHORT32 FLAG! */
#include "../include/config.h"
#include "../include/int32type.h"
#include <assert.h>
/* For doing byte reversal on 32-bit words */

void byterevn(int32type w[],int n)
{
  register int32type old;
  register int32type newv;
  int j;
  1?((void )0) : __assert_fail("sizeof(int32type) == 4","byterevn.c",16,"void byterevn(int *, int)");
  for (j = 0; j < n; j++) {
    old = w[j];
    newv = ((old >> 24) & 0x000000ff);
    newv |= ((old >> 8) & 0x0000ff00);
    newv |= ((old << 8) & 0x00ff0000);
    newv |= ((old << 24) & 0xff000000);
    w[j] = newv;
  }
/* byterevn */
}
