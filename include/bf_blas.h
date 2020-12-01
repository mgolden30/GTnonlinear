#ifndef BF_BLAS
#define BF_BLAS

#include "bf_work_mem.h"

void bf_blas_swap( int n, bf *x, int incx, bf *y, int incy );
void bf_blas_scal( int n, bf alpha, bf *x, int incx );
void bf_blas_copy( int n, bf *x, int incx, bf *y, int incy );
void bf_blas_axpy(int n, bf alpha, bf *x, int incx, bf *y, int incy );
void bf_blas_dot( bf result, int n, bf *x, int incx, bf *y, int incy );

void qr_decomposition( bf *qt, bf *a, int m, int n, int lda );
void bf_print_matrix(bf *matrix, int m, int n, int lda);
#endif
