/******** spectrum2.c *************/
/* MIMD version 6 */
/* Spectrum for Kogut-Susskind pointlike hadrons, wall source 
* MIMD version 6 with SG additions
* updated 2/21/95, by SG.
* modified 7/96 DT, do both Naik and conventional spectrum, wall
*	sources only. 
* 7/29/96 DT fixed names of rho and rho2 propagators.
*   rho is VT (rho-b_1), rho2 is PV (rho-a_1)
*   also pi and pi2 -> pi_ps_prop, pi_sc_prop
* modified 8/96 DT changed source and source_inc
* modified 5/97 DT for improved actions
* modified 2/99 DT for new output format
*
* This version does arbitrary number of wall sources 
* This version DOES NOT fix the gauge -- you should do that before 
*	calling it
*
*/
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"
*/

#include "generic_ks_includes.h"

#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif

int spectrum2( float vmass, field_offset temp1, field_offset temp2 ){ 
  /* return the C.G. iteration number */
  double *pi_ps_prop,*pi_sc_prop,*rho_pv_prop,*rho_vt_prop,*barprop;
  float vmass_x2;
  register complex cc;
  float finalrsq;
  int isrc;
  register int i,x,y,z,t,icol,cgn;
  register int t_source,t_off;

  vmass_x2 = 2.*vmass;
  cgn=0;
  pi_ps_prop = (double *)malloc(nt*sizeof(double) );   /* "pi" */
  pi_sc_prop = (double *)malloc(nt*sizeof(double) );   /* "pi2" */
  rho_vt_prop = (double *)malloc(nt*sizeof(double) );  /* "rho" */
  rho_pv_prop = (double *)malloc(nt*sizeof(double) );  /* "rho2" */
  barprop = (double *)malloc(nt*sizeof(double) );
  for( t=0; t<nt; t++){
    pi_ps_prop[t] = pi_sc_prop[t] = rho_pv_prop[t] = rho_vt_prop[t]=0.0;
    barprop[t] = 0.0;
  }

  for(t_source=source_start, isrc=0; t_source<2*nt && isrc < n_sources; 
	++isrc, t_source += source_inc ) {
	if(this_node==0)printf("spectrum(): source time = %d\n",t_source);

	  for(icol=0; icol<3; icol++) {
		      
	      /* initialize temp1 and temp2 */
	      clear_latvec( temp1, EVENANDODD);
	      clear_latvec( temp2, EVENANDODD);

	      if(t_source < nt){  /* wall source */
	         for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
		     if( node_number(x,y,z,t_source) != mynode() )continue;
		     i=node_index(x,y,z,t_source);
		     ((su3_vector *)(F_PT(&lattice[i],temp1)))->c[icol].real
		       = -1.0;
	 	 }
	      }
		
	      /* do a C.G. (source in temp1, result in temp2) */
	      if(t_source%2 == 0) {
		cgn += ks_congrad( temp1, temp2, vmass,
		    niter,rsqprop,EVEN,&finalrsq);
	          /* Multiply by -Madjoint */
	          dslash( temp2, F_OFFSET(ttt), ODD);
	          scalar_mult_latvec( temp2, -vmass_x2, F_OFFSET(ttt),
			EVEN);
/**copy_latvec( temp1, F_OFFSET(g_rand), EVENANDODD );
checkmul();**/
/**check_invert( F_OFFSET(ttt), temp1 );**/
	      }
	      else {
		cgn += ks_congrad( temp1, temp2, vmass,
		    niter,rsqprop,ODD,&finalrsq);
	          /* Multiply by -Madjoint */
	          dslash( temp2, F_OFFSET(ttt), EVEN);
	          scalar_mult_latvec( temp2, -vmass_x2, F_OFFSET(ttt),
			ODD);
	      }
	      
	      /* fill the hadron matrix */
	      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);
	    } /* end loop on icol */
	  
	  /* measure the meson propagator */
	  for(t=0; t<nt; t++){
	      /* define the time value offset t from t_source */
	      t_off = (t+t_source)%nt;
	      
	      for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)
		for(icol=0;icol<3;icol++) {
		    if( node_number(x,y,z,t_off) != mynode() )continue;
		    i=node_index(x,y,z,t_off);
		    cc = su3_dot( &lattice[i].propmat[icol],
				 &lattice[i].propmat[icol] );
		    
		    pi_ps_prop[t] += cc.real;
		    
		    if( (x+y)%2==0)rho_pv_prop[t] += cc.real;
		    else	   rho_pv_prop[t] -= cc.real;
		    if( (y+z)%2==0)rho_pv_prop[t] += cc.real;
		    else	   rho_pv_prop[t] -= cc.real;
		    if( (z+x)%2==0)rho_pv_prop[t] += cc.real;
		    else	   rho_pv_prop[t] -= cc.real;
		    
		    if( x%2==0)rho_vt_prop[t] += cc.real;
		    else       rho_vt_prop[t] -= cc.real;
		    if( y%2==0)rho_vt_prop[t] += cc.real;
		    else       rho_vt_prop[t] -= cc.real;
		    if( z%2==0)rho_vt_prop[t] += cc.real;
		    else       rho_vt_prop[t] -= cc.real;
		    
		    if( (x+y+z)%2==0)pi_sc_prop[t] += cc.real;
		    else	     pi_sc_prop[t] -= cc.real;
		    
		  } /* icol */
	      
	    } /* nt-loop */
	  
	  /* measure the baryon propagator */
	  for(t=0; t<nt; t++) {
	      /* define the time value offset t from t_source */
	      t_off = (t+t_source)%nt;
	      
	      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
		  if( node_number(x,y,z,t_off) != mynode() )continue;
		  i=node_index(x,y,z,t_off);
		  cc = det_su3( (su3_matrix *)lattice[i].propmat );
			/* propmat[3] is su3_vectors, needs typecast */
		      
		  /* must get sign right.  This looks to see if we have
			wrapped around the lattice.  "t" is the distance
			from the source to the measurement, so we are
			trying to find out if t_source+t is greater than
			or equal to nt.  the "-tsource/nt" is in there
			so that it will work correctly with tsource=nt.  */
		  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
		     barprop[t] += cc.real;
		  else  /* change sign because antiperiodic b.c.  sink point
			should really be in a copy of the lattice */
		     barprop[t] -= cc.real;
	      }
	    } /* nt-loop */
  
    } /* end loop on t_source */

  /* dump the propagators */
  g_vecdoublesum( pi_ps_prop, nt );
  g_vecdoublesum( rho_pv_prop, nt );
  g_vecdoublesum( rho_vt_prop, nt );
  g_vecdoublesum( pi_sc_prop, nt );
  g_vecdoublesum( barprop, nt );
  if( this_node==0 ){
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",vmass,vmass);
    printf("SOURCE: CORNER\n");
    printf("SINKS: PION_PS PION_SC RHO_VT RHO_PV\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0 %e 0.0 %e 0.0 %e 0.0\n",t,
      pi_ps_prop[t]/n_sources, pi_sc_prop[t]/n_sources,
      rho_vt_prop[t]/n_sources, rho_pv_prop[t]/n_sources);
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e  %e\n",vmass,vmass,vmass);
    printf("SOURCE: CORNER\n");
    printf("SINKS: NUCLEON\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0\n",t,
      barprop[t]/n_sources);
    printf("ENDPROP\n");

    fflush(stdout);
  }
  free( pi_ps_prop ); free( pi_sc_prop );
  free( rho_pv_prop ); free( rho_vt_prop );
  free( barprop );
  return(cgn);
} /* spectrum */

