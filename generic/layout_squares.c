/******** layout_squares.c *********/
/* MIMD version 6 */
/* UNTESTED!! C. DeTar 7/26/99 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version puts two dimensional slices on nodes.  The slices are
   taken in the two shortest directions, so that the areas of the
   interfaces between slices will be as small as possible.

   The remaining two dimensions are divided up into rectangles.  We
   require that the product of the two longest dimensions is divisible
   by the number of nodes, a power of two.

   3/29/00 EVENFIRST is the rule now. CD.
*/

/*
   setup_layout() does any initial setup.  When it is called the
     lattice dimensions nx,ny,nz and nt have been set.
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

int dirs[4];	/* list of directions, longest first */
int dims[4];	/* lattice dimensions, in order X,Y,Z,T */
int nsites_per;	/* number of sites per plane */
int xsquaresize,ysquaresize;	/* dimensions of rectangle */
int nxsquares,nysquares;	/* number of rectangles in each direction */

void setup_layout(){
register int i,j,k;
    if(mynode()==0){
	printf("LAYOUT = 2d-squares, options = ");
	printf("EVENFIRST,");
	printf("\n");
    }

    /* sort directions with longest ones first */
    dirs[0]=XUP; dirs[1]=YUP; dirs[2]=ZUP; dirs[3]=TUP;
    dims[XUP]=nx; dims[YUP]=ny; dims[ZUP]=nz; dims[TUP]=nt;
    for(i=3;i>0;i--)for(j=0;j<i;j++){
	if(dims[dirs[j]]<dims[dirs[j+1]]){
	    k=dirs[j]; dirs[j]=dirs[j+1]; dirs[j+1]=k;
	}
    }
    nsites_per = dims[dirs[ZUP]]*dims[dirs[TUP]];

    /* Figure out dimensions of rectangle */
    i=1;	/* current number of squares */
    xsquaresize = dims[dirs[XUP]];
    ysquaresize = dims[dirs[YUP]];
    nxsquares = nysquares = 1;
    while(i<numnodes()){
	if( dims[dirs[XUP]]/nxsquares >= dims[dirs[YUP]]/nysquares ){
	    if( (dims[dirs[XUP]]/nxsquares)%2 == 0){ /* decrease x size */
		i*=2; xsquaresize /= 2; nxsquares *= 2;
	    }
	    else if( (dims[dirs[YUP]]/nysquares)%2 == 0){ /* dec.. y size */
		i*=2; ysquaresize /= 2; nysquares *= 2;
	    }
	    else {
		if(mynode()==0)printf(
	       "LAYOUT: Can't lay out this lattice, not enough factors of 2\n");
		terminate(1);
	    }
	}
	else {
	    if( (dims[dirs[YUP]]/nysquares)%2 == 0){ /* decrease y size */
		i*=2; ysquaresize /= 2; nysquares *= 2;
	    }
	    else if( (dims[dirs[XUP]]/nxsquares)%2 == 0){ /* dec. x size */
		i*=2; xsquaresize /= 2; nxsquares *= 2;
	    }
	    else {
		if(mynode()==0)printf(
	       "LAYOUT: Can't lay out this lattice, not enough factors of 2\n");
		terminate(1);
	    }
	}
    }

    if( mynode()==0)if( dims[dirs[ZUP]]%2 != 0 || dims[dirs[TUP]]%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
    sites_on_node = nsites_per*(xsquaresize*ysquaresize);
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(int x,int y,int z,int t) {
register int i,xr,yr;
int coords[4];
    coords[XUP]=x; coords[YUP]=y; coords[ZUP]=z; coords[TUP]=t;
    xr = coords[dirs[XUP]]/xsquaresize;
    yr = coords[dirs[YUP]]/ysquaresize;
    i = xr + nxsquares*yr;
    return( i );
}

int node_index(int x,int y,int z,int t) {
register int i,xr,yr;
int coords[4];
    coords[XUP]=x; coords[YUP]=y; coords[ZUP]=z; coords[TUP]=t;
    xr = coords[dirs[XUP]]%xsquaresize;
    yr = coords[dirs[YUP]]%ysquaresize;
    i = nsites_per*(xr + xsquaresize*yr) + coords[dirs[ZUP]] +
	dims[dirs[ZUP]]*coords[dirs[TUP]];
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + nsites_per*xsquaresize*ysquaresize)/2 );
    }
}

int num_sites(int node) {
    return( nsites_per*xsquaresize*ysquaresize );
}
