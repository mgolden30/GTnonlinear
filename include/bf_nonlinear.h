#ifndef BFNONLIN
#define BFNONLIN

#include "bf_blas.h"

//Takes R^n -> R^m
typedef struct bf_nonlinear{
    int m;
    int n;
    void (*function)(bf *out, bf *in, bf *numeric_params, void *params );
    bf    *numeric_params;
    void  *params;
}bf_nonlinear;

//use this macro to set y = F(x)
#define BF_NONLINEAR_EVAL(y,F,x)   ( (*(F.function))(y,x,(F.numeric_params),(F.params)) )

void newton_raphson( bf_nonlinear f, bf *state, int max_iterations, bf hookstep, bf threshold, bf tolerance);

#endif
