#include "bf_nonlinear.h"
#include "bf_dynamics.h"

void planar_3body_td(bf *out, bf *in);
void periodic_orbit_objective_function( bf *out, bf *in, bf *np, void *p );

//Symmetries
void fix_rotation( bf *state, bf tolerance );
void rescale_Hamiltonian( bf *state );


int main( int argc, char *argv[] ){
    bfp prec = 300;
    int m = 8;
    int n = 9;

    bf error;
    bf_init( error, prec );

    bf threshold, tolerance, hookstep;
    bf_init( threshold, prec );
    bf_init( tolerance, prec );
    bf_init( hookstep,  prec );
    bf_set_d( tolerance, 1e-50 );
    bf_set_d( error, 1e-15);
 
    bf_nonlinear f;
    f.m = m;
    f.n = n;
    f.function = &periodic_orbit_objective_function;
    f.numeric_params = NULL;
    f.params = &error;

    bf state[n];
    bf_inits( n, state, prec );

    double state_d[] = {
     1.7193284422e+00,
-5.7337446177e-01,
 1.3995529123e-01,
 8.3754768588e-02,
-2.6178729333e-01,
 7.6679641899e-01,
-1.2127842743e-01,
-1.1860925952e-01,
 1.7082388113e+01
    };
    for(int i=0; i<n; i++){
        bf_set_d( state[i], state_d[i] );
    }

    bf_print_vector(n, state);
   
    printf("Starting Newton.\n"); 
   
    bf_set_d( hookstep, 0.001 );
    bf_set_d( threshold, 1e-3 );
    newton_raphson( f, state, 1000, hookstep, threshold, tolerance);
   
   /* 
    bf_set_d( hookstep, 1.0 );
    bf_set_d( threshold, 1e-9 );
    newton_raphson( f, state, 10, hookstep, threshold, tolerance);
   */

    bf final_state[8];
    bf_inits( 8, final_state, prec );

    rescale_Hamiltonian( state );
    int output = 0; //For output to a file
    rk45( final_state, state, state[8], error, &planar_3body_td, 8, output );

    bf_print_vector( n, state );

    return 0;
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
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(in[0]);
    check_work_mem( &work, 9, prec);
    
    bf *copy_of_in = &work.ptr[0];
    bf *tolerance  = &work.ptr[8];
    bf_set_d( *tolerance, 1e-35 );

    bf_blas_copy( 8, in, 1, copy_of_in, 1 );

    int dimension = 8;
    //in[8] is the time to evolve until
    bf *error    = (bf *) p;
    int oodt     = 500;
    int output = 1; //No output during newton
    rk45( out, in, in[8], *error, &planar_3body_td, dimension, output );
    
    //fix_rotation( out, *tolerance );
    //fix_rotation( copy_of_in, *tolerance );
    
    for(int i=0; i<8; i++){
        bf_sub( out[i], out[i], copy_of_in[i] );
    }
}

/* Takes in a 8D state and rotates it so that r1 is on the x-axis.
 *
 */
void fix_rotation( bf *state, bf tolerance ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(state[0]);
    check_work_mem( &work, 4, prec);
 
    bf *q = work.ptr;

    qr_decomposition( q, state, 2, 1, 1, tolerance );

    bf_blas_mv(2, 2, q, 2, &state[2], 1, &state[2], 1);
    bf_blas_mv(2, 2, q, 2, &state[4], 1, &state[4], 1);
    bf_blas_mv(2, 2, q, 2, &state[6], 1, &state[6], 1);
}

//Takes in a 9-d state and rescales it so that H=-1
void rescale_Hamiltonian( bf *state ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(state[0]);
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
    
    //Now fill r_i
    bf_set(r1[0], state[0]);
    bf_set(r1[1], state[1]);
    bf_set(r2[0], state[2]);
    bf_set(r2[1], state[3]);
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
    
    //Compute the Hamiltonian resuing old memory
    bf *temp        = &r1[0];
    bf *hamiltonian = &r1[1];
    bf_set_ui( *hamiltonian, 0);

    //Add all the kinetic energy
    for(int i=4; i<8; i++){
        bf_blas_dot( *temp, 1, &state[i], 1, &state[i], 1 );
	bf_div_ui(*temp, *temp, 2);
	bf_add( *hamiltonian, *hamiltonian, *temp );
    }
    
    //Now the potential energy
    bf_ui_div(*d12, 1, *d12);
    bf_ui_div(*d13, 1, *d13);
    bf_ui_div(*d23, 1, *d23);

    bf_sub( *hamiltonian, *hamiltonian, *d12 );
    bf_sub( *hamiltonian, *hamiltonian, *d13 );
    bf_sub( *hamiltonian, *hamiltonian, *d23 );
    
    //Now compute the scaling parameter and store it in hamiltonian

    bf *lambda = hamiltonian;
    //Take the negative hamiltonian so the sqrt is well-defined
    bf_neg( *lambda, *lambda );
    bf_sqrt( *lambda, *lambda );
    bf_ui_div( *lambda, 1, *lambda );

    // r -> lambda^-2 r
    for( int i=0; i<4; i++){
        bf_div( state[i], state[i], *lambda);
        bf_div( state[i], state[i], *lambda);
    }
    // p -> lambda * p
    for( int i=4; i<8; i++){
        bf_mul( state[i], state[i], *lambda);
    }
    //time -> lambda^-3 * time 
    bf_div( state[8], state[8], *lambda );
    bf_div( state[8], state[8], *lambda );
    bf_div( state[8], state[8], *lambda );
    
    //All done! Your dynamics are preserved and H=-1
}
