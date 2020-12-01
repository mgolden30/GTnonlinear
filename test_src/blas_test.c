#include "bf_blas.h"

int main( int argc, char *argv[] ){
    int n=4;
    bf x[n];
    bf y[n];
    bfp prec = 100;
    bf_inits(n, x, prec);
    bf_inits(n, y, prec);
    for(int i=0; i<n; i++){
        bf_set_ui(x[i], i);
        bf_set_ui(y[i], 4*i+3);
    }
 
    printf("x = \n");
    bf_print_vector(n, x);
    printf("y = \n");
    bf_print_vector(n, y);

    bf dot;
    bf_init( dot, prec );
    bf_blas_dot(dot, n, x, 1, y, 1);
    printf("x^T y = ");
    bf_print(dot);
    printf("\n");


    //Now test QR decomposition
    printf("Testing QR decomposition\n"); 

    n=3;
    //let x from above be our matrix
    bf qt[n*n];
    bf r[n*n];
    bf_inits(n*n, qt, prec);
    bf_inits(n*n, r, prec);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            bf_set_ui( r[i*n+j], i*n+j+1);
	}
    }


    bf_print_matrix( r, n, n, n);
    
	    qr_decomposition( qt, r, n, n, n );

    bf_print_matrix( qt, n, n, n);
    bf_print_matrix( r, n, n, n);
    
    return 0;
}
