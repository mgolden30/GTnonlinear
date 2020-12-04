#include "bf_nonlinear.h"

void rk4( bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n );
void make_timeseries( bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n );
void rk23( bf *out, bf *in, bf end_time, bf error, void (*td)(bf *out, bf *in), int n );
void rk45( bf *out, bf *in, bf end_time, bf error, void (*td)(bf *out, bf *in), int n );
void planar_3body_td(bf *out, bf *in);
void periodic_orbit_objective_function( bf *out, bf *in, bf *np, void *p );




int main( int argc, char *argv[] ){
    bfp prec = 300;
    int m = 8;
    int n = 9;

    bf error;
    bf_init( error, prec );
    bf_set_d(error, 1e-10);

    bf_nonlinear f;
    f.m = m;
    f.n = n;
    f.function = &periodic_orbit_objective_function;
    f.numeric_params = NULL;
    f.params = &error;

    bf state[n];
    bf_inits( n, state, prec );

    //PO1
/*    double state_d[] = {9.3213883162e-01, 
 4.1184551339e-01, 
-1.7388921759e-01, 
-4.5095085400e-02,
-3.5597373723e-01,
 7.8417112947e-01,
-3.1838720613e-01,
 6.4203522194e-01,
 1.4135554594e+01};*/


    //PO3
/*    double state_d[] = {4.991812150770e-01,
-1.195770031071e-01,
-5.505054841634e-01,
 1.699918044280e+00,
 7.487473342824e-02,
-7.233755902280e-01,
 1.262530333574e-01,
 2.782645088666e-01,
 4.959828632869e+00}; */


    //PO2
    double state_d[] = {8.444489449743e-02,
 2.153332620005e-02,
 8.369080946320e-01,
 7.341808005541e-01,
 7.535202858588e-01,
-8.938040908186e-01,
-6.552806166016e-01,
 7.348976526846e-01,
 7.327423125946e+00};

    for(int i=0; i<n; i++){
        bf_set_d( state[i], state_d[i] );
    }

    bf_print_vector(n, state);

    bf threshold, tolerance, hookstep;
    bf_init( threshold, prec );
    bf_init( tolerance, prec );
    bf_init( hookstep,  prec );
    bf_set_d( threshold, 1e-6 );
    bf_set_d( tolerance, 1e-50 );
    bf_set_d( hookstep, 0.1 );
   
    printf("Starting Newton.\n"); 
    newton_raphson( f, state, 100, hookstep, threshold, tolerance);

    bf final_state[8];
    bf_inits( 8, final_state, prec );
    make_timeseries( final_state, state, state[8], 500, &planar_3body_td, 8);

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


void make_timeseries( bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n ){
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

    FILE *orbit = fopen("orbit.dat", "w");
    bf_set_ui( *time, 0 );
    while( bf_cmp(*time, end_time) < 0 ){
	for(int i=0; i<8; i++){
	    mpfr_fprintf( orbit, "% .10Rf\t", out[i] );
        }
	fprintf(orbit, "\n");

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

    fclose(orbit);
}


//Adaptive Runge-Kutta using the midpoint method/Ralston's third order method
void rk23( bf *out, bf *in, bf end_time, bf error, void (*td)(bf *out, bf *in), int n ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(in[0]);
    check_work_mem( &work, 6*n+3, prec);
    /* k1 - n
     * k2 - n
     * k3 - n
     * temp - n
     * 
     * d2 - n //change predicted by RK2
     * d3 - n //change predicted by RK3
     *
     * alpha - 1
     * time - 1
     * dt   - 1
     */
    bf *k1    = &work.ptr[0];
    bf *k2    = &work.ptr[n];
    bf *k3    = &work.ptr[2*n];
    bf *temp  = &work.ptr[3*n];

    bf *d2    = &work.ptr[4*n];
    bf *d3    = &work.ptr[5*n];

    bf *alpha = &work.ptr[6*n];
    bf *time  = &work.ptr[6*n+1];
    bf *dt    = &work.ptr[6*n+2];

    //start by copying the current state to out
    bf_blas_copy( n, in, 1, out, 1 );

    bf_set_ui( *time, 0 );
    bf_set_ui( *dt, 1);
    while( bf_cmp(*time, end_time) < 0 ){
	//See if time + dt > end_time
        bf_add( temp[0], *time, *dt);
       	if( bf_cmp( temp[0], end_time ) > 0 ){
            bf_sub( *dt, end_time, *time );
	}

compute_time_deriv:
        td(k1, out);
        bf_div_ui(*alpha, *dt, 2);
        bf_blas_copy( n, out, 1, temp, 1 );
	bf_blas_axpy( n, *alpha, k1, 1, temp, 1);

	td(k2, temp);
        bf_blas_copy( n, out, 1, temp, 1 );
        bf_mul_ui(*alpha, *dt,    3);
        bf_div_ui(*alpha, *alpha, 4);
	bf_blas_axpy( n, *alpha, k2, 1, temp, 1);

	td(k3, temp);

        //Now compute predictions from both integrators
	//alpha is repurposed here to be a rational number * dt
        for(int i=0; i<n; i++){
            bf_set_ui( d2[i], 0 );
            bf_set_ui( d3[i], 0 );
	}
	
	bf_blas_axpy( n, *dt, k2, 1, d2, 1);

	bf_div_ui( *alpha, *dt,    9 );
	bf_mul_ui( *alpha, *alpha, 2 );
	bf_blas_axpy( n, *alpha, k1, 1, d3, 1);
        bf_div_ui( *alpha, *dt, 3 );
	bf_blas_axpy( n, *alpha, k2, 1, d3, 1);
        bf_div_ui( *alpha, *dt,    9 );
	bf_mul_ui( *alpha, *alpha, 4 );
	bf_blas_axpy( n, *alpha, k3, 1, d3, 1);
        
        //Now use alpha to compute the distance between these two increments.
	for(int i=0; i<n; i++){
            bf_sub( d2[i], d3[i], d2[i] );
	}
	bf_blas_dot( *alpha, n, d2, 1, d2, 1);
	bf_sqrt( *alpha, *alpha );
     
       /*	
 	printf("error estimate is ");
        bf_print( *alpha );
	printf(" with a timestep of ");
        bf_print( *dt );
	printf("\n");
*/


	if( bf_cmp(*alpha, error) > 0){
            //timestep is too big apparently.
	    //cut it in half and try again
	    bf_div_ui(*dt, *dt, 2);
	    goto compute_time_deriv;
	}

        //If you have made it this far, you have a sufficiently small step
	//adjust your state by d3, which is the step predicted by RK3
	for(int i=0; i<n; i++){
            bf_add(out[i], out[i], d3[i]);
	}

        //Increment time!
	bf_add(*time, *time, *dt);

	//Increase timestep by a percent.
	bf_mul_d(*dt, *dt, 1.01);
    }
}



// The Runge–Kutta–Fehlberg method. It is a fifth order integrator with an embedded fourth order integrator so that 
// error can be estimated. An adaptive step-size can be used.
void rk45( bf *out, bf *in, bf end_time, bf error, void (*td)(bf *out, bf *in), int n ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(in[0]);
    check_work_mem( &work, 9*n+3, prec);
    /* k1 - n
     * k2 - n
     * k3 - n
     * k4 - n
     * k5 - n
     * k6 - n
     *
     * temp - n
     * 
     * d4 - n //change predicted by 4th order method
     * d5 - n //change predicted by 5th order method
     *
     * alpha - 1
     * time - 1
     * dt   - 1
     */

    bf *k1    = &work.ptr[0*n];
    bf *k2    = &work.ptr[1*n];
    bf *k3    = &work.ptr[2*n];
    bf *k4    = &work.ptr[3*n];
    bf *k5    = &work.ptr[4*n];
    bf *k6    = &work.ptr[5*n];

    bf *temp  = &work.ptr[6*n];

    bf *d4    = &work.ptr[7*n];
    bf *d5    = &work.ptr[8*n];

    bf *alpha = &work.ptr[9*n];
    bf *time  = &work.ptr[9*n+1];
    bf *dt    = &work.ptr[9*n+2];

    //start by copying the current state to out
    bf_blas_copy( n, in, 1, out, 1 );

    bf_set_ui( *time, 0 );
    bf_set_ui( *dt, 1);
    while( bf_cmp(*time, end_time) < 0 ){
	//See if time + dt > end_time
        bf_add( temp[0], *time, *dt);
       	if( bf_cmp( temp[0], end_time ) > 0 ){
            bf_sub( *dt, end_time, *time );
	}

compute_time_deriv:
        td(k1, out);

        bf_blas_copy( n, out, 1, temp, 1 );
        bf_div_ui(*alpha, *dt, 4);
	bf_blas_axpy( n, *alpha, k1, 1, temp, 1);

	td(k2, temp);

        bf_blas_copy( n, out, 1, temp, 1 );
        bf_mul_ui(*alpha, *dt,     3);
        bf_div_ui(*alpha, *alpha, 32);
	bf_blas_axpy( n,  *alpha, k1, 1, temp, 1);
        bf_mul_ui(*alpha, *dt,     9);
        bf_div_ui(*alpha, *alpha, 32);
	bf_blas_axpy( n,  *alpha, k2, 1, temp, 1);

	td(k3, temp);

        bf_blas_copy( n, out, 1, temp, 1 );
        bf_mul_ui(*alpha, *dt,    1932);
        bf_div_ui(*alpha, *alpha, 2197);
	bf_blas_axpy( n,  *alpha, k1, 1, temp, 1);
        bf_mul_si(*alpha, *dt,   -7200);
        bf_div_ui(*alpha, *alpha, 2197);
	bf_blas_axpy( n,  *alpha, k2, 1, temp, 1);
        bf_mul_ui(*alpha, *dt,    7296);
        bf_div_ui(*alpha, *alpha, 2197);
	bf_blas_axpy( n,  *alpha, k3, 1, temp, 1);

        td(k4, temp);

        bf_blas_copy( n, out, 1, temp, 1 );
        bf_mul_ui(*alpha, *dt,    439);
        bf_div_ui(*alpha, *alpha, 216);
	bf_blas_axpy( n,  *alpha, k1, 1, temp, 1);
        bf_mul_si(*alpha, *dt,   -8);
	bf_blas_axpy( n,  *alpha, k2, 1, temp, 1);
        bf_mul_ui(*alpha, *dt,    3680);
        bf_div_ui(*alpha, *alpha, 513);
	bf_blas_axpy( n,  *alpha, k3, 1, temp, 1);
        bf_mul_si(*alpha, *dt,    -845);
        bf_div_ui(*alpha, *alpha, 4104);
	bf_blas_axpy( n,  *alpha, k4, 1, temp, 1);

        td(k5, temp);

        bf_blas_copy( n, out, 1, temp, 1 );
        bf_mul_si(*alpha, *dt,    -8);
        bf_div_ui(*alpha, *alpha, 27);
	bf_blas_axpy( n,  *alpha, k1, 1, temp, 1);
        bf_mul_ui(*alpha, *dt,   2);
	bf_blas_axpy( n,  *alpha, k2, 1, temp, 1);
        bf_mul_si(*alpha, *dt,    -3544);
        bf_div_ui(*alpha, *alpha, 2565);
	bf_blas_axpy( n,  *alpha, k3, 1, temp, 1);
        bf_mul_ui(*alpha, *dt,    1859);
        bf_div_ui(*alpha, *alpha, 4104);
	bf_blas_axpy( n,  *alpha, k4, 1, temp, 1);
        bf_mul_si(*alpha, *dt,    -11);
	bf_div_ui(*alpha, *alpha, 40);
	bf_blas_axpy( n,  *alpha, k5, 1, temp, 1);

        td(k6, temp);

        //Now compute predictions from both integrators
	//alpha is repurposed here to be a rational number * dt
        for(int i=0; i<n; i++){
            bf_set_ui( d4[i], 0 );
            bf_set_ui( d5[i], 0 );
	}
	
	//First compute d5
	bf_mul_ui(*alpha, *dt,    16);
        bf_div_ui(*alpha, *alpha, 135);
	bf_blas_axpy( n,  *alpha, k1, 1, d5, 1);
	bf_mul_ui(*alpha, *dt,    6656);
        bf_div_ui(*alpha, *alpha, 12825);
	bf_blas_axpy( n,  *alpha, k3, 1, d5, 1);
	bf_mul_ui(*alpha, *dt,    28561);
        bf_div_ui(*alpha, *alpha, 56430);
	bf_blas_axpy( n,  *alpha, k4, 1, d5, 1);
	bf_mul_si(*alpha, *dt,    -9);
        bf_div_ui(*alpha, *alpha, 50);
	bf_blas_axpy( n,  *alpha, k5, 1, d5, 1);
	bf_mul_ui(*alpha, *dt,    2);
        bf_div_ui(*alpha, *alpha, 55);
	bf_blas_axpy( n,  *alpha, k6, 1, d5, 1);

        //Now compute d4
        bf_mul_ui(*alpha, *dt,    25);
        bf_div_ui(*alpha, *alpha, 216);
	bf_blas_axpy( n,  *alpha, k1, 1, d4, 1);
        bf_mul_ui(*alpha, *dt,    1408);
        bf_div_ui(*alpha, *alpha, 2565);
	bf_blas_axpy( n,  *alpha, k3, 1, d4, 1);
        bf_mul_ui(*alpha, *dt,    2197);
        bf_div_ui(*alpha, *alpha, 4104);
	bf_blas_axpy( n,  *alpha, k4, 1, d4, 1);
        bf_mul_si(*alpha, *dt,    -1);
        bf_div_ui(*alpha, *alpha, 5);
	bf_blas_axpy( n,  *alpha, k5, 1, d4, 1);

        //Now use alpha to compute the distance between these two increments.
	for(int i=0; i<n; i++){
            bf_sub( d4[i], d5[i], d4[i] );
	}

	bf_blas_dot( *alpha, n, d4, 1, d4, 1);
	bf_sqrt( *alpha, *alpha );
/*
  	printf("error estimate is ");
        bf_print( *alpha );
	printf(" with a timestep of ");
        bf_print( *dt );
	printf("\n");
*/      

        if( bf_cmp(*alpha, error) > 0){
            //timestep is too big apparently.
	    //cut it in half and try again
	    bf_div_ui(*dt, *dt, 2);
	    goto compute_time_deriv;
	}

        //If you have made it this far, you have a sufficiently small step
	//adjust your state by d3, which is the step predicted by RK3
	for(int i=0; i<n; i++){
            bf_add(out[i], out[i], d5[i]);
	}

        //Increment time!
	bf_add(*time, *time, *dt);

	//Increase timestep by a percent.
	bf_mul_d(*dt, *dt, 1.01);
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
    int dimension = 8;
    //in[8] is the time to evolve until
    bf *error = (bf *) p;
    int oodt = 500;
    rk45( out, in, in[8], *error, &planar_3body_td, dimension );
    //rk23( out, in, in[8], *error, &planar_3body_td, dimension );
    //rk4( out, in, in[8], oodt, &planar_3body_td, dimension );
    for(int i=0; i<8; i++){
        bf_sub( out[i], out[i], in[i] );
    }
}

