#include "bf_nonlinear.h"


/* Performs Newton-Raphson iteration on f to find a root r: f(r) = 0
 *
 * IN:
 * bf_nonlinear f     - a nonlinear function you want to find a root of
 * bf *state          - the state you start with, hopefully near a root. 
 * int max_iterations - so you won't wander state-space forever
 * bf hookstep        - a big_float between 0 and 1 that weights your Newton step.
 *                      a value closer to zero will converge more slowly, but prevents overstepping.
 * bf threshold       - exit when ||f(state)|| < threshold
 * bf tolerance       - the size of a value we take to be negligible. Needed for numerical linear algebra
 *                      if a value is smaller than tolerance, it is assumed to be from numerical error.
 */
void newton_raphson( bf_nonlinear f, bf *state, int max_iterations, bf hookstep, bf threshold, bf tolerance){
    /* Need the following working memory
     * jacobian         - m*n   
     * jacobian_t       - m*n   tranpose of jacobian
     * perturbed vector - n
     * image            - m
     * perturb_image    - m
     * q                - n*n
     * step             - n
     * h                - 1
     */

    int m = f.m;
    int n = f.n;

    BF_WORKING_MEMORY(work);
    bfp prec = bf_get_prec(state[0]);
    check_work_mem( &work, 2*m*n + n+ 2*m + n*n + n + 1 , prec);
    bf *jacobian   = &work.ptr[0*m*n];
    bf *jacobian_t = &work.ptr[1*m*n];
    bf *perturb    = &work.ptr[2*m*n];
    bf *image      = &work.ptr[2*m*n+n];
    bf *pert_image = &work.ptr[2*m*n+n+m];
    bf *q          = &work.ptr[2*m*n+n+2*m];
    bf *step       = &work.ptr[2*m*n+n+2*m+n*n];
    bf *h          = &work.ptr[2*m*n+n+2*m+n*n+n];


    for( int iteration=0; iteration < max_iterations; iteration++ ){
        //Make h.
	//I will choose for h to decrease with the size of f(state)
	BF_NONLINEAR_EVAL(image, f, state);
	bf_blas_dot( *h, m, image, 1, image, 1);
	bf_sqrt(*h,*h);

	//take this oppurtunity to check if the state has converged.
        printf("Iteration %03d: error is ", iteration);
        bf_print(*h);
        printf("\n");
	
	if( bf_cmp( *h, threshold ) <= 0){
	    printf("State converged!\n");
            return;
        }
	bf_div_ui(*h, *h, 1000000000);
        //Now h = ||f(state)||*1E-9


        //Fill out the jacobian.
	for( int i=0; i<n; i++){
            bf_blas_copy( n, state, 1, perturb, 1);
	    bf_add(perturb[i], perturb[i], *h);
            BF_NONLINEAR_EVAL( pert_image, f, perturb);
            for(int j=0; j<m; j++){
                bf_sub( jacobian[i + j*n], pert_image[j], image[j] );
		bf_div( jacobian[i + j*n], jacobian[i + j*n], *h );
	    }
	}


	//Now that jacobian is done, take its transpose.
        bf_matrix_transpose( jacobian, m, n, jacobian_t);

        //Now plug it into QR decomposition to effectively do LQ decomposition
        qr_decomposition( q, jacobian_t, n, m, m, tolerance );
	//printf("Rotated Jacobian L^T: where J = LQ\n");
        //bf_print_matrix( jacobian_t, n, m, m);       

	//Take transpose of q
	for(int i=0; i<n; i++){
            for(int j=i+1; j<n; j++){
                bf_swap( q[i*n+j], q[j*n+i] );
	    }  
	}

	//Now invert jacobian_t transpose. It is a left triangular matrix.
        bf_div( step[0], image[0], jacobian_t[0]); //assume this first step can be done without trouble.
	for(int i=1; i<m; i++){
	    //use h since it isn't needed anymore.
	    bf_blas_dot(*h, i, &jacobian_t[i], m, step, 1);

	    bf_sub( *h, image[i], *h);
	    if( bf_cmp( jacobian_t[i*m+i], threshold ) <= 0){
                bf_set_ui( step[i], 0 );
                continue;
	    }
	    bf_div( step[i], *h, jacobian_t[i*m+i] );
        }
        //Set remaining spots to zero
	for(int i=m; i<n; i++){
            bf_set_ui( step[i], 0 );
	}

	/*
        printf("I have now found a rotated step vector.\n");
        
	bf l[m*n];
	bf_inits( m*n, l, prec);

        bf_matrix_transpose( jacobian_t, n, m, l);
        printf("Here is L:\n");
	bf_print_matrix( l, m, n, n);

        printf("Here is the image:\n");
	bf_print_vector( m, image );
        
        printf("Here is l*s':\n");
        bf_blas_mv( m,  n, l, n, step, 1, image, 1);
	bf_print_vector( m, image );
        
        printf("Now let us check that the jacobian can be faithfully reproduced.\n");
	printf("J = \n");
	bf_print_matrix( jacobian, m, n, n);
	printf("L * q = \n");
        bf_blas_mm( m, n, n, l, n, q, n, jacobian, n);
	bf_print_matrix( jacobian, m, n, n);
        */

	/*
        printf("Before taking the step, let us verify our step actually will resolve the error.\n");
	bf_print_vector( m, image );
	*/
	//and just multiply it by q!
        bf_blas_mv( n, n, q, n, step, 1, step, 1);
	
	/*
	printf("\nHere is Jacobian*step.\n");
        bf_blas_mv( m, n, jacobian, n, step, 1, image, 1);
	bf_print_vector( m, image );*/

        //Now take the step.
        for( int i=0; i<n; i++ ){
            bf_mul( step[i], step[i], hookstep );
	    bf_sub( state[i], state[i], step[i] );
	}
    }
}
