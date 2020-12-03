#include "bf_nonlinear.h"

void rk4( bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n );
void planar_3body_td(bf *out, bf *in);
void periodic_orbit_objective_function( bf *out, bf *in, bf *np, void *p );




int main( int argc, char *argv[] ){
    bfp prec = 200;
    int m = 8;
    int n = 9;

    bf_nonlinear f;
    f.m = m;
    f.n = n;
    f.function = &periodic_orbit_objective_function;
    f.numeric_params = NULL;
    f.params = NULL;

    bf state[n];
    bf_inits( n, state, prec );
    double state_d[] = { .9305, 0.4161, -0.1719, -0.0499, -0.3584, 0.7829, -0.3170, 0.6432, 14.1377 };
    for(int i=0; i<n; i++){
        bf_set_d( state[i], state_d[i] );
    }

    bf_print_vector(n, state);

    bf threshold, tolerance, hookstep;
    bf_init( threshold, prec );
    bf_init( tolerance, prec );
    bf_init( hookstep,  prec );
    bf_set_d( threshold, 1e-4 );
    bf_set_d( tolerance, 1e-10 );
    bf_set_d( hookstep,  0.2 );
   
    printf("Starting Newton.\n"); 
    newton_raphson( f, state, 100, hookstep, threshold, tolerance);

    bf_print_vector( n, state );

    return 0;
}





/* Runge-Kutta 4th order
 *
 * IN:
 * bf *out - the result of time integration
 * bf *in  - the initial state
 * bf end_time
 * unsigned int oodt - one over dt
 * void (*td) time derivative function.
 * int n - dimension of the space.
 */
void rk4( bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(in[0]);
    check_work_mem( &work, 5*n+3, prec);
    /* k1 - n
     * k2 - n
     * k3 - n
     * k4 - n
     * temp - n
     *
     * alpha - 1
     * time - 1
     * dt   - 1
     */
    bf *k1    = &work.ptr[0];
    bf *k2    = &work.ptr[n];
    bf *k3    = &work.ptr[2*n];
    bf *k4    = &work.ptr[3*n];
    bf *temp  = &work.ptr[4*n];
    bf *alpha = &work.ptr[5*n];
    bf *time  = &work.ptr[5*n+1];
    bf *dt    = &work.ptr[5*n+2];

    //start by copying the current state to out
    bf_blas_copy( n, in, 1, out, 1 );

    bf_set_ui( *time, 0 );
    while( bf_cmp(*time, end_time) < 0 ){
        bf_set_ui(*dt, 1);
	bf_div_ui(*dt, *dt, oodt);
	//By default set dt to 1/oodt
	
	bf_add( temp[0], *time, *dt); //See if time + dt > end_time
       	if( bf_cmp( temp[0], end_time ) > 0 ){
            bf_sub( *dt, end_time, *time );
	}

        td(k1, out);
        bf_div_ui(*alpha, *dt, 2);
        bf_blas_copy( n, out, 1, temp, 1 );
	bf_blas_axpy( n, *alpha, k1, 1, temp, 1);

	td(k2, temp);
        bf_blas_copy( n, out, 1, temp, 1 );
	bf_blas_axpy( n, *alpha, k2, 1, temp, 1);

	td(k3, temp);
        bf_blas_copy( n, out, 1, temp, 1 );
	bf_blas_axpy( n, *dt, k3, 1, temp, 1);

	td(k4, temp);

	//Now add these to out to increment the state.
        bf_div_ui( *alpha, *dt, 6 );
	bf_blas_axpy( n, *alpha, k1, 1, out, 1);
	bf_blas_axpy( n, *alpha, k4, 1, out, 1);
        bf_div_ui( *alpha, *dt, 3 );
	bf_blas_axpy( n, *alpha, k2, 1, out, 1);
	bf_blas_axpy( n, *alpha, k3, 1, out, 1);

        //Increment time!
	bf_add(*time, *time, *dt);
    }
}


void planar_3body_td(bf *out, bf *in){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(in[0]);
    check_work_mem( &work, 15, prec);
    /* r1  - 2
     * r2  - 2
     * r3  - 2
     * r12 - 2
     * r13 - 2
     * r23 - 2
     * d12 - 1
     * d13 - 1
     * d23 - 1
     */
    bf *r1  = &work.ptr[0];
    bf *r2  = &work.ptr[2];
    bf *r3  = &work.ptr[4];
    bf *r12 = &work.ptr[6];
    bf *r13 = &work.ptr[8];
    bf *r23 = &work.ptr[10];
    bf *d12 = &work.ptr[12];
    bf *d13 = &work.ptr[13];
    bf *d23 = &work.ptr[14];
    
    //First use the momenta
    bf_set( out[0], in[4] );
    bf_set( out[1], in[5] );
    bf_set( out[2], in[6] );
    bf_set( out[3], in[7] );

    //Now fill r_i
    bf_set(r1[0], in[0]);
    bf_set(r1[1], in[1]);
    bf_set(r2[0], in[2]);
    bf_set(r2[1], in[3]);
    bf_add(r3[0], r1[0], r2[0]);
    bf_add(r3[1], r1[1], r2[1]);
    bf_neg(r3[0], r3[0]);
    bf_neg(r3[1], r3[1]);

    //Take the differences
    bf_sub( r12[0], r1[0], r2[0] );
    bf_sub( r12[1], r1[1], r2[1] );
    bf_sub( r13[0], r1[0], r3[0] );
    bf_sub( r13[1], r1[1], r3[1] );
    bf_sub( r23[0], r2[0], r3[0] );
    bf_sub( r23[1], r2[1], r3[1] );

    //Compute the magnitudes of these vectors.
    bf_blas_dot( *d12, 2, r12, 1, r12, 1);
    bf_blas_dot( *d13, 2, r13, 1, r13, 1);
    bf_blas_dot( *d23, 2, r23, 1, r23, 1);
    bf_sqrt( *d12, *d12 );
    bf_sqrt( *d13, *d13 );
    bf_sqrt( *d23, *d23 );

    //Divide r_ij by the appropriate power of magnitudes
    for(int i=0; i<2; i++){
        bf_div( r12[i], r12[i], *d12 );
        bf_div( r12[i], r12[i], *d12 );
        bf_div( r12[i], r12[i], *d12 );

        bf_div( r13[i], r13[i], *d13 );
        bf_div( r13[i], r13[i], *d13 );
        bf_div( r13[i], r13[i], *d13 );

        bf_div( r23[i], r23[i], *d23 );
        bf_div( r23[i], r23[i], *d23 );
        bf_div( r23[i], r23[i], *d23 );
    }

    for(int i=4; i<8; i++){
        bf_set_ui( out[i], 0 );
    } 

    for(int i=0; i<2; i++){
        bf_sub( out[4+i], out[4+i], r12[i] );
        bf_sub( out[4+i], out[4+i], r13[i] );
        bf_add( out[6+i], out[6+i], r12[i] );
        bf_sub( out[6+i], out[6+i], r23[i] );
    }
}



void periodic_orbit_objective_function( bf *out, bf *in, bf *np, void *p ){
    rk4( out, in, in[8], 200, &planar_3body_td, 8);
    for(int i=0; i<8; i++){
        bf_sub( out[i], out[i], in[i] );
    }
}

