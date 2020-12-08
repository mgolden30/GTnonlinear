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
    check_work_mem( &work, m, prec );
    
    for(int i=0; i<m; i++){
        bf_blas_dot( work.ptr[i], n, &a[i*lda], 1, x, incx );
    }

    //Now swap the values of y and working memory.
    //We do this so that the possibility x=y is allowed.
    for(int i=0; i<m; i++){
        bf_set( y[i*incy], work.ptr[i] );
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
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(matrix[0]);
    check_work_mem( &work, m*n, prec);
 
    for( int i=0; i<m; i++ ){
        for( int j=0; j<n; j++ ){
            bf_set( work.ptr[i+j*m], matrix[j+i*n] );
	}
    }
    bf_blas_copy(m*n, work.ptr, 1, matrix_t, 1 );
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
	//also check for NAN errors. I am seeing these when I use QR to diagonalize a matrix
	if( bf_cmp( *alpha, tolerance ) <= 0  ){
            continue;
	}

	//Find mu
        bf_sub( *mu, *alpha, a_vec[0] );
	bf_mul( *mu, *mu, *alpha );
	bf_mul_ui( *mu, *mu, 2 );
	bf_sqrt( *mu, *mu );

        //Check that mu is well-defined currently
        //If mu is too small, NANs will appear and you will be a sad camper
	if(  bf_cmp( *mu, tolerance ) <= 0  ){
            continue;
	}

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

    //Explicitly set small values of a to zero no.
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
	    //store absolute value in alpha
	    bf_abs( *alpha, a[i*lda+j] );
            if( bf_cmp(*alpha, tolerance) <= 0 ){
               bf_set_ui( a[i*lda+j], 0);
	    }
	}
    }
}


/* Diagonalize a symmetric, positive definite (SPD) matrix.
 * Uses the QR algorithm
 */
void bf_diagonalize_spd( bf *u, bf *a, int m, int lda, bf tolerance){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(a[0]);
    check_work_mem( &work, m*m, prec);
    
    bf *rotation = &work.ptr[0];

    //Set u to the identity
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
	    unsigned int d = (i==j) ? 1 : 0; 
            bf_set_ui( u[i*m+j], d);
        }
    }

    for(int i=0; i<100; i++){
        qr_decomposition( rotation, a, m, m, lda, tolerance );
	printf("rotation derived by qr is \n");
        bf_print_matrix( rotation, m, m, m);

        bf_blas_mm( m, m, m, rotation, m, u, m, u, m);

	bf_matrix_transpose( rotation, m, m, rotation );
        bf_blas_mm( m, m, m, a, lda, rotation, m, a, lda);

	printf("iteration %d:\n", i);
	printf("D\n");
	bf_print_matrix( a, m, m, m );
	printf("U\n");
	bf_print_matrix( u, m, m, m );
    }
    //When you are done, we need to take the transpose of u
    bf_matrix_transpose( u, m, m, u );
}


/* Finds the singular value decomposition of a matrix
 */
void bf_svd( bf *u, bf *vt, bf *sigma, bf *a, int m, int n, int lda, bf tolerance){
    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(a[0]);
    check_work_mem( &work, m*m+n*n+m*n, prec);
    
    bf *aat = &work.ptr[0];
    bf *ata = &work.ptr[m*m];
    bf *at  = &work.ptr[m*m+n*n];

    //Set at to the transpose of a
    //This does not respect lda for now
    bf_matrix_transpose( a, m, n, at);    
    bf_blas_mm( m, n, m, a,  lda, at, m,   aat, m);
    bf_blas_mm( n, m, n, at, m,   a,  lda, ata, n);
    
    bf_diagonalize_spd( u,  aat, m, m, tolerance );
    bf_diagonalize_spd( vt, ata, m, m, tolerance );

    //Now u is correct, but vt contains v.
    //Sigma = U^T A V, so do the V multiplication now
    bf_blas_mm( m, n, n, a, lda, vt, n, sigma, n);
    //Take transpose of vt and store transpose of u in work since we don't need any of it anymore
    bf *ut = &work.ptr[0];
    bf_matrix_transpose(u,  m, m, ut);
    bf_matrix_transpose(vt, n, n, vt);
    
    //finish finding Sigma!
    bf_blas_mm( m, m, n, u, m, sigma, n, sigma, n);


    printf("SVD accomplished.\n");
    printf("U = \n");
    bf_print_matrix( u, m, m, m );
    printf("Sigma = \n");
    bf_print_matrix( sigma, m, n, n );
    printf("V^T = \n");
    bf_print_matrix( vt, n, n, n );
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
