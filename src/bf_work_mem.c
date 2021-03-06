#include "big_float.h"
#include "stdlib.h"

typedef struct bf_work_mem{
    bf *ptr;
    bfp prec;
    int size; //Size of the array ptr points to
}bf_work_mem;

#define BF_WORKING_MEMORY(name)  static bf_work_mem name = {NULL, 0, 0}

//Check to make sure there is enough memory allocated
void check_work_mem( bf_work_mem *w, int n, bfp p ){
    if( w->size < n ){
        //Not enough memory is currently allocated
	if( w->ptr != NULL ) free(w->ptr);
	w->ptr  = (bf*) malloc(n*sizeof(bf));
	w->size = n;
	w->prec = p;
        bf_inits(n, w->ptr, p);
    }
}
