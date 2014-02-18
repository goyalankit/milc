/******** layout_hyper_tstretch.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */
/* This version divides the lattice by factors of prime numbers in any of the
   four directions.  It prefers to divide the longest dimensions,
   which mimimizes the area of the surfaces.  Similarly, it prefers
   to divide dimensions which have already been divided, thus not
   introducing more off-node directions.
	S. Gottlieb, May 18, 1999
	The code will start trying to divide with the largest prime factor
	and then work its way down to 2.  The current maximum prime is 53.
	The array of primes on line 46 may be extended if necessary.
   This requires that the lattice volume be divisible by the number
   of nodes.  Each dimension must be divisible by a suitable factor
   such that the product of the four factors is the number of nodes.
   DT 1/7/03  Allow last hypercubes in time direction to have one less
   time slice than others.  For example, put a 24^3x64 lattice on 120
   processors by dividing time direction as 4*13+12.  This is NOT implemented
   in the most general way.
   3/29/00 EVENFIRST is the rule now. CD.
   12/10/00 Fixed so k = MAXPRIMES-1 DT
*/
/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
     This routine sets the global variables "sites_on_node",
     "even_sites_on_node" and "odd_sites_on_node".
   num_sites(node) returns the number of sites on a node
   node_number(x,y,z,t) returns the node number on which a site lives.
   node_index(x,y,z,t) returns the index of the site on the node - ie the
     site is lattice[node_index(x,y,z,t)].
   These routines will change as we change our minds about how to distribute
     sites among the nodes.  Hopefully the setup routines will work for any
     consistent choices. (ie node_index should return a different value for
     each site on the node.)
*/
#include "generic_includes.h"
/* dimensions of hypercubes */
int squaresize[4UL];
/* number of hypercubes in each direction */
int nsquares[4UL];
int prime[] = {(2), (3), (5), (7), (11), (13), (17), (19), (23), (29), (31), (37), (41), (43), (47), (53)};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )

void setup_layout()
{
  register int i;
  register int j;
  register int k;
  register int dir;
  if (mynode() == 0) {
    printf("LAYOUT = Hypercubes, options = ");
    printf("EVENFIRST,");
    printf("\n");
  }
/* Figure out dimensions of rectangle */
  if (((((nx * ny) * nz) * nt) % numnodes()) == 0) {
    squaresize[0] = nx;
    squaresize[1] = ny;
    squaresize[2] = nz;
    squaresize[3] = nt;
  }
  else if (((((nx * ny) * nz) * (nt + 1)) % numnodes()) == 0) {
    squaresize[0] = nx;
    squaresize[1] = ny;
    squaresize[2] = nz;
    squaresize[3] = (nt + 1);
  }
  nsquares[0] = (nsquares[1] = (nsquares[2] = (nsquares[3] = 1)));
/* current number of hypercubes */
  i = 1;
  while(i < numnodes()){
/* figure out which prime to divide by starting with largest */
    k = (sizeof(prime) / sizeof(int ) - 1);
    while((((numnodes() / i) % prime[k]) != 0) && (k > 0))
      --k;
/* figure out which direction to divide */
/* find largest even dimension of h-cubes */
    for (((j = 1) , (dir = 0)); dir <= 3; dir++) 
      if ((squaresize[dir] > j) && ((squaresize[dir] % prime[k]) == 0)) 
        j = squaresize[dir];
{
/* if one direction with largest dimension has already been
	   divided, divide it again.  Otherwise divide first direction
	   with largest dimension. */
      for (dir = 0; dir <= 3; dir++) 
        if ((squaresize[dir] == j) && (nsquares[dir] > 1)) 
          break; 
    }
    if (dir > 3) {
      for (dir = 0; dir <= 3; dir++) 
        if (squaresize[dir] == j) 
          break; 
    }
/* This can fail if I run out of prime factors in the dimensions */
    if (dir > 3) {
      if (mynode() == 0) 
        printf("LAYOUT: Can\'t lay out this lattice, not enough factors of %d\n",prime[k]);
      terminate(1);
    }
/* do the surgery */
    i *= prime[k];
    squaresize[dir] /= prime[k];
    nsquares[dir] *= prime[k];
  }
  sites_on_node = (((squaresize[0] * squaresize[1]) * squaresize[2]) * squaresize[3]);
  if (!(((((nx * ny) * nz) * nt) % numnodes()) == 0) && (((((nx * ny) * nz) * (nt + 1)) % numnodes()) == 0)) {
/* stretched t direction by one */
    if (mynode() == 0) 
      printf("SOME NODES HAVE FEWER SITES\n");
    if (mynode() >= (((nsquares[0] * nsquares[1]) * nsquares[2]) * (nsquares[3] - 1))) {
      sites_on_node = (((squaresize[0] * squaresize[1]) * squaresize[2]) * (squaresize[3] - 1));
    }
  }
/* Need even number of sites per hypercube */
  if (mynode() == 0) 
    if ((sites_on_node % 2) != 0) {
      printf("SORRY, CAN\'T LAY OUT THIS LATTICE\n");
      terminate(0);
    }
  if (mynode() == 0) 
    printf("ON EACH NODE %d x %d x %d x %d\n",squaresize[0],squaresize[1],squaresize[2],squaresize[3]);
  if ((mynode() == 0) && ((sites_on_node % 2) != 0)) 
    printf("WATCH OUT FOR EVEN/ODD SITES ON NODE BUG!!!\n");
  even_sites_on_node = (odd_sites_on_node = (sites_on_node / 2));
}

int node_number(int x,int y,int z,int t)
{
  register int i;
  x /= squaresize[0];
  y /= squaresize[1];
  z /= squaresize[2];
  t /= squaresize[3];
  i = (x + (nsquares[0] * (y + (nsquares[1] * (z + (nsquares[2] * t))))));
  return i;
}

int node_index(int x,int y,int z,int t)
{
  register int i;
  register int xr;
  register int yr;
  register int zr;
  register int tr;
  xr = (x % squaresize[0]);
  yr = (y % squaresize[1]);
  zr = (z % squaresize[2]);
  tr = (t % squaresize[3]);
  i = (xr + (squaresize[0] * (yr + (squaresize[1] * (zr + (squaresize[2] * tr))))));
/* even site */
  if (((((x + y) + z) + t) % 2) == 0) {
    return i / 2;
  }
  else {
    return (i + sites_on_node) / 2;
  }
}

int num_sites(int node)
{
  return sites_on_node;
}
