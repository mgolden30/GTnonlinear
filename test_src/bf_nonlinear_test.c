#include "bf_nonlinear.h"

//R^3 -> R^2
void func( bf *out, bf *in, bf *np, void *p){
    bf_mul( out[0], in[0], in[1] );
    bf_sub_ui( out[0], out[0], 1 );
    bf_mul( out[1], in[0], in[2] );
    bf_sub_ui( out[1], out[1], 2 );
}

int main( int argc, char *argv[] ){
    bfp prec = 200;
    int m = 2;
    int n = 3;

    bf_nonlinear f;
    f.m = m;
    f.n = n;
    f.function = &func;
    f.numeric_params = NULL;
    f.params = NULL;

    bf state[n];
    bf_inits( n, state, prec );
    for(int i=0; i<n; i++){
        bf_set_ui( state[i], i+1 );
    }

    bf threshold, tolerance, hookstep;
    bf_init( threshold, prec );
    bf_init( tolerance, prec );
    bf_init( hookstep,  prec );
    bf_set_d( threshold, 1e-30 );
    bf_set_d( tolerance, 1e-40 );
    bf_set_d( hookstep,  0.2 );
    
    newton_raphson( f, state, 400, hookstep, threshold, tolerance);

    bf_print_vector( n, state );

    return 0;
}
