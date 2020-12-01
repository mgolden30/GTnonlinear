#include <stdio.h>
#include <stdlib.h>

typedef struct gmp_function{
    int in_dim;
    int out_dim;
    void (*function)(mpf_t *out, mpf_t *in, mpf_t *numeric_params, void *params, gmp_work work);
    gmp_work work;
    mpf_t *numeric_params;
    void  *params;
}gmp_function;

//use this macro to set y = F(x)
#define GMP_FN_EVAL(y,F,x)   ( (*(F.function))(y,x,(F.numeric_params),(F.params),(F.work)) )



void identity(mpf_t *out, mpf_t *in, mpf_t *numeric_params, void *params){
    mpf_set( *out, *in );
}

int main(int argc, char *argv[] ){
    gmp_function id;
    id.in_dim  = 1;
    id.out_dim = 1;
    id.function = &identity;
    id.numeric_params = NULL;
    id.params = NULL;

    mpf_t x,y;
    const int prec = 100;
    mpf_init2(x, prec);
    mpf_init2(y, prec);

    mpf_set_ui(x, 17);

    GMP_FN_EVAL(&y, id, &x);

    mpf_out_str(stdout, 10, 20, y);

    mpf_clears(x,y);
    return 0;
}
