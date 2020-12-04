#include "bf_blas.h"

int main( int argc, char *argv[] ){
    int n=4;
    bf x[n];
    bf y[n];
    bfp prec = 200;
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
    printf("Testing QR decomposition\n\n"); 

    int m=4; n=5;
    //let x from above be our matrix
    bf qt[m*m];
    bf r[m*n];
    bf_inits(m*m, qt, prec);
    bf_inits(m*n, r,  prec);
    for(int i=0; i<m*n; i++){    
        bf_set_ui( r[i], i*(i+1));
    }

    bf tolerance;
    bf_init( tolerance, prec );
    bf_set_d( tolerance, 1e-20 );

    bf_print_matrix( r, m, n, n);
    
    qr_decomposition( qt, r, m, n, n, tolerance );

    bf_print_matrix( qt, m, m, m);
    bf_print_matrix( r,  m, n, n);

    //reset r
    for(int i=0; i<m*n; i++){    
        bf_set_ui( r[i], i*(i+1));
    }
    
    //multiply by qt to check that the triangular matrix is found
    bf_print_matrix( r,  m, n, n);
    bf_blas_mm( m,  m, n, qt, m, r, n, r, n);
    bf_print_matrix( r,  m, n, n);
 
    return 0;
}
