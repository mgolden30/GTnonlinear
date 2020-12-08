#ifndef BF_DYNAMICS
#define BF_DYNAMICS

#include "bf_blas.h"

//Standard Runge-Kutta
void rk4(  bf *out, bf *in, bf end_time, unsigned int oodt, void (*td)(bf *out, bf *in), int n );

//Adaptive Runge-Kutta
void rk23( bf *out, bf *in, bf end_time,          bf error, void (*td)(bf *out, bf *in), int n );
void rk45( bf *out, bf *in, bf end_time,          bf error, void (*td)(bf *out, bf *in), int n, int output );

#endif
