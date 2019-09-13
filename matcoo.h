#ifndef SPARSEMAT_MATCOO_H
#define SPARSEMAT_MATCOO_H
#include "ll_float.h"
#include <stdio.h>

typedef struct matcoo_ matcoo;
struct matcoo_ {
    ll_float* coords_val;
    ll_float* coords_i;
    ll_float* coords_j;
    int dimX;
    int dimY;
};

matcoo* matcoo_new(const float* data, int dimX, int dimY);

#endif //SPARSEMAT_MATCOO_H