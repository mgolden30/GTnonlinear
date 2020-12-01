#ifndef BFWORK
#define BFWORK

#include "big_float.h"
#include "stdlib.h"

typedef struct bf_work_mem{
    bf *ptr;
    bfp prec;
    int size; //Size of the array ptr points to
}bf_work_mem;

#define BF_WORKING_MEMORY(name)  static bf_work_mem name = {NULL, 0, 0}

void check_work_mem( bf_work_mem *w, int n, bfp p );

#endif
