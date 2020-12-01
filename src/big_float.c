#include "big_float.h"

void bf_init( bf x,  bfp prec ){
    mpfr_init2( x, prec );
}

void bf_inits( int n, bf *x, bfp prec){
    for(int i=0; i<n; i++){
        bf_init(x[i], prec);
    }
}

void bf_clear( bf x ){
    mpfr_clear(x);
}

void bf_clears( int n, bf *x ){
    for(int i=0; i<n; i++){
        bf_clear(x[i]);
    }
}

bfp bf_get_prec( bf x){
    return mpfr_get_prec(x);
}

//rop <- op
void bf_set( bf rop, bf op ){
    mpfr_set( rop, op, MPFR_RNDN );
}

//rop <- op
void bf_set_ui( bf rop, unsigned int op ){
    mpfr_set_ui( rop, op, MPFR_RNDN );
}

//rop <- op
void bf_set_d( bf rop, double op ){
    mpfr_set_d( rop, op, MPFR_RNDN );
}

void bf_swap( bf x, bf y){
    mpfr_swap(x, y);
}


void bf_neg(bf rop, bf op){
    mpfr_neg( rop, op, MPFR_RNDN );
}

//rop <- op1 + op2
void bf_add(bf rop, bf op1, bf op2){
    mpfr_add( rop, op1, op2, MPFR_RNDN );
}

//rop <- op1 - op2
void bf_sub(bf rop, bf op1, bf op2){
    mpfr_sub( rop, op1, op2, MPFR_RNDN );
}

//rop <- op1 * op2
void bf_mul(bf rop, bf op1, bf op2){
    mpfr_mul( rop, op1, op2, MPFR_RNDN );
}

//rop <- op1 * op2
void bf_mul_ui(bf rop, bf op1, unsigned int op2){
    mpfr_mul_ui( rop, op1, op2, MPFR_RNDN );
}

//rop <- op1 / op2
void bf_div(bf rop, bf op1, bf op2){
    mpfr_div( rop, op1, op2, MPFR_RNDN );
}

//rop <- op1 / op2
void bf_ui_div(bf rop, unsigned int op1, bf op2){
    mpfr_ui_div( rop, op1, op2, MPFR_RNDN );
}


//rop <- op1 / op2
void bf_div_ui(bf rop, bf op1, unsigned int op2){
    mpfr_div_ui( rop, op1, op2, MPFR_RNDN );
}

void bf_sqrt(bf rop, bf op){
    mpfr_sqrt( rop, op, MPFR_RNDN );
}

void bf_print(bf x){
    mpfr_out_str( stdout, 10, 6, x, MPFR_RNDN );
}

void bf_print_vector(int n, bf *x){
    for( int i=0; i<n; i++){
        bf_print(x[i]);
	printf("\n");
    }
    printf("\n");
}

