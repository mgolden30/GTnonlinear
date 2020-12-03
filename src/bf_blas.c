/* This file implements many useful linear algebra concepts for big_floats
 * we will store matrices in row major order.
 */

#include "bf_blas.h"

void bf_blas_swap( int n, bf *x, int incx, bf *y, int incy ){
    for(int i=0; i<n; i++){
        bf_swap( x[i*incx], y[i*incy] );
    }
}

void bf_blas_scal( int n, bf alpha, bf *x, int incx ){
    for(int i=0; i<n; i++){
        bf_mul(x[i*incx], alpha, x[i*incx]);
    }
}

//y <- x
void bf_blas_copy( int n, bf *x, int incx, bf *y, int incy ){
    for(int i=0; i<n; i++){
        bf_set( y[i*incy], x[i*incx] );
    }
}

//y <- y + ax
void bf_blas_axpy(int n, bf alpha, bf *x, int incx, bf *y, int incy){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(x[0]);
    check_work_mem( &work, 1, prec );

    for( int i=0; i<n; i++ ){
        bf_mul(work.ptr[0], alpha, x[i*incx]);
	bf_add( y[i*incy], work.ptr[0], y[i*incy] );
    }
}

//result = x^T y
void bf_blas_dot( bf result, int n, bf *x, int incx, bf *y, int incy ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(x[0]);
    check_work_mem( &work, 1, prec );

    bf_set_ui( result, 0 );
    for( int i=0; i<n; i++){
        bf_mul( work.ptr[0], x[i*incx], y[i*incy] );
	bf_add( result, result, work.ptr[0] );
    }
}


/*
 * Level 2 BLAS
 */

//I made this because gemv is a pain in the ass to call.
//Just normal matrix-vector multiplication
//y<-a*x
void bf_blas_mv(int m, int n, bf *a, int lda, bf *x, int incx, bf *y, int incy){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(a[0]);
    check_work_mem( &work, n, prec );
    
    for(int i=0; i<m; i++){
        bf_blas_dot( work.ptr[i], n, &a[i*lda], 1, x, incx );
    }

    //Now swap the values of y and working memory.
    //We do this so that the possibility x=y is allowed.
    for(int i=0; i<m; i++){
        bf_swap( work.ptr[i], y[i*incy] );
    }
}



/*
 * Level 3 BLAS
 */
//I made this because gemm is a pain in the ass to call.
//Just normal matrix-matrix multiplication
//c<-a*b
void bf_blas_mm( int m, int n, int k, bf *a, int lda, bf *b, int ldb, bf *c, int ldc){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(a[0]);
    check_work_mem( &work, m*k, prec );
    
    for(int i=0; i<m; i++){
        for(int j=0; j<k; j++){
            bf_blas_dot( work.ptr[i*k+j], n, &a[i*lda], 1, &b[j], ldb );
	}
    }

    //Now swap the values with working memory.
    //We do this so that the possibility c=a or c=b is allowed.
    for(int i=0; i<m; i++){
        for(int j=0; j<k; j++){
            bf_swap( work.ptr[i*k+j], c[i*ldc+j] );
        }
    }
}




/* Take the transpose of a matrix.
 */
void bf_matrix_transpose( bf *matrix, int m, int n, bf *matrix_t ){
    for( int i=0; i<m; i++ ){
        for( int j=0; j<n; j++ ){
            bf_set( matrix_t[i+j*m], matrix[j+i*n] );
	}
    }
}





/* Here are some more complicated linear algebra functions
 */


/* Perofrms a QR decomposition of a matrix.
 * A = QR, where A is a general m-by-n matrix, Q is a m-by-m orthogonal matrix, and R is an upper triangular m-by-n matrix.
 * 
 * IN:
 * bf *qt - a pointer to Q^T. Must be initialized, but not set to anything in particular.
 * bf *a  - matrix A you want to do the QR decomposition for. This is destroyed in the process.
 * int m,n -size of bf *a
 * int lda - this is a Lapack-style variable in case you want to perform QR on a submatrix. Usually lda = n
 * bf tolerance - a big_float that determines what is numerical error.
 *
 * OUT:
 * bf *qt - now has the transpose of q in it.
 * bf *a  - now has the upper-triangular matrix.
 */
void qr_decomposition( bf *qt, bf *a, int m, int n, int lda, bf tolerance ){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(a[0]);
    check_work_mem( &work, m+2, prec);
   
    //See "Unitary triangularization of a nonzymmetric matrix" by Alston Householder to see how this is done
    //I have tried to keep a similar notation to him.
    bf *alpha = &work.ptr[0];
    bf *mu    = &work.ptr[1];
    bf *u_vec = &work.ptr[2];

    //Start off by setting qt to the identity
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
	    unsigned int d = (i==j) ? 1 : 0; 
            bf_set_ui( qt[i*m+j], d);
	}
    }

    //Decide how many householder transformations we need
    int num_hh = m>n ? n : m-1;
    for(int i=0; i<num_hh; i++){
        bf *a_vec        = &a[i*lda + i]; //The ith column, starting at the ith element
	int a_vec_length = m-i;

	//Find alpha
	bf_blas_dot( *alpha, a_vec_length, a_vec, lda, a_vec, lda);
	bf_sqrt( *alpha, *alpha );
        
	//Check if the norm of this vector is less than your tolerance
	//If it is, then this vector is effectively already zero. Skip to the next column
	if( bf_cmp( *alpha, tolerance ) <= 0 ){
            continue;
	}

	//Find mu
        bf_sub( *mu, *alpha, a_vec[0] );
	bf_mul( *mu, *mu, *alpha );
	bf_mul_ui( *mu, *mu, 2 );
	bf_sqrt( *mu, *mu );
	//turn a_vec into u_vec
       	bf_blas_copy( a_vec_length, a_vec, lda, u_vec, 1);
        bf_sub( u_vec[0], u_vec[0], *alpha );
        //turn mu into 1/mu and scale the vector.
	bf_ui_div( *mu, 1, *mu);
	bf_blas_scal( a_vec_length, *mu, u_vec, 1 );

	for(int j=i; j<n; j++){
            a_vec = &a[i*lda + j];
	    //store dot product in alpha since it isn't needed again.
	    bf_blas_dot( *alpha, a_vec_length, a_vec, lda, u_vec, 1 );
	    bf_neg(*alpha, *alpha);
	    bf_mul_ui( *alpha, *alpha, 2);
            bf_blas_axpy( a_vec_length, *alpha, u_vec, 1, a_vec, lda);
	}

        //Do the same routine to qt
	for(int j=0; j<m; j++){
            bf *q_vec = &qt[i*m + j];
	    //store dot product in alpha since it isn't needed again.
	    bf_blas_dot( *alpha, a_vec_length, q_vec, m, u_vec, 1 );
	    bf_neg(*alpha, *alpha);
	    bf_mul_ui( *alpha, *alpha, 2);
            bf_blas_axpy( a_vec_length, *alpha, u_vec, 1, q_vec, m);
	}
    }
}

void bf_print_matrix(bf *matrix, int m, int n, int lda){
    for(int j=0; j<m; j++){
        for(int i=0; i<n; i++){
            bf_print(matrix[j*n+i]);
	    printf("  ");
	}
	printf("\n");
    }
    printf("\n");
}
