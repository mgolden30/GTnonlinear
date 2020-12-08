#ifndef BF_BLAS
#define BF_BLAS

#include "bf_work_mem.h"

void bf_blas_swap( int n, bf *x, int incx, bf *y, int incy );
void bf_blas_scal( int n, bf alpha, bf *x, int incx );
void bf_blas_copy( int n, bf *x, int incx, bf *y, int incy );
void bf_blas_axpy(int n, bf alpha, bf *x, int incx, bf *y, int incy );
void bf_blas_dot( bf result, int n, bf *x, int incx, bf *y, int incy );

void bf_matrix_transpose( bf *matrix, int m, int n, bf *matrix_t );

void bf_blas_mv(int m, int n, bf *a, int lda, bf *x, int incx, bf *y, int incy);

void bf_blas_mm( int m, int n, int k, bf *a, int lda, bf *b, int ldb, bf *c, int ldc);



void qr_decomposition( bf *qt, bf *a, int m, int n, int lda, bf tolerance );
void bf_diagonalize_spd( bf *u, bf *a, int m, int lda, bf tolerance);
void bf_svd( bf *u, bf *vt, bf *sigma, bf *a, int m, int n, int lda, bf tolerance);

void bf_print_matrix(bf *matrix, int m, int n, int lda);
#endif
