#ifndef BIG_FLOAT
#define BIG_FLOAT

#include <stdio.h>
#include <mpfr.h>

typedef mpfr_t     big_float;
typedef big_float  bf;

typedef mpfr_prec_t         big_float_precision;
typedef big_float_precision bfp;


void bf_init( bf x,  bfp prec );
void bf_inits( int n, bf *x, bfp prec);
void bf_clear( bf x );
void bf_clears( int n, bf *x );
bfp bf_get_prec( bf x);
void bf_set( bf rop, bf op );
void bf_set_ui( bf rop, unsigned int op );
void bf_set_d( bf rop, double op );
void bf_swap( bf x, bf y);
void bf_neg(bf rop, bf op);
int  bf_cmp(bf op1, bf op2);
void bf_add(bf rop, bf op1, bf op2);
void bf_sub(bf rop, bf op1, bf op2);
void bf_sub_ui(bf rop, bf op1, unsigned int op2);
void bf_mul(bf rop, bf op1, bf op2);
void bf_mul_ui(bf rop, bf op1, unsigned int op2);
void bf_div(bf rop, bf op1, bf op2);
void bf_ui_div(bf rop, unsigned int op1, bf op2);
void bf_div_ui(bf rop, bf op1, unsigned int op2);
void bf_sqrt(bf rop, bf op);

void bf_print(bf x);
void bf_print_vector(int n, bf *x);

#endif
