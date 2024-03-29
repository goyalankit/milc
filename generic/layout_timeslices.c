/******** layout_timeslices.c *********/
/* MIMD version 6 */
/* ROUTINES WHICH DETERMINE THE DISTRIBUTION OF SITES ON NODES */

/* This version puts entire timeslices on nodes. It requires that nt,
   the time extent, is a multiple of the number of nodes used.
   We hope this speeds up spatial FFT's

   3/29/00 EVENFIRST is the rule now. CD.
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

int nslices;		/* number of timeslices per node */

void setup_layout(){
    if(mynode()==0){
	printf("LAYOUT = Timeslices, options = ");
	printf("EVENFIRST,");
	printf("\n");
    }

    if( nt%numnodes() !=0 ){
	if(mynode()==0)printf(
	    "LAYOUT: Can't lay out this lattice: nt not multiple of nummodes\n");
	    terminate(1);
	}

    nslices = nt / numnodes();
    sites_on_node = nx*ny*nz*nslices;
    /* Need even number of sites per hypercube */
    if( mynode()==0)if( sites_on_node%2 != 0){
	printf("SORRY, CAN'T LAY OUT THIS LATTICE\n");
	terminate(0);
    }
    even_sites_on_node = odd_sites_on_node = sites_on_node/2;
}

int node_number(int x, int y, int z, int t) {
    t /= nslices;
    return( t );
}

int node_index(int x, int y, int z, int t) {
register int i, tr;
    tr = t % nslices;
    i = x + nx*(y + ny*(z + nz*tr));
    if( (x+y+z+t)%2==0 ){	/* even site */
	return( i/2 );
    }
    else {
	return( (i + sites_on_node)/2 );
    }
}

int num_sites(int node) {
    return( sites_on_node );
}
