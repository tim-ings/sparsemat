#ifndef SPARSEMAT_MATCOO_H
#define SPARSEMAT_MATCOO_H
#include "ll_float.h"
#include <stdio.h>

typedef struct matcoo_ matcoo;
struct matcoo_ {
    ll_float* vals;
    ll_float* is;
    ll_float* js;
    int dimX;
    int dimY;
};

matcoo* matcoo_new(const float* data, int dimX, int dimY);
void matcoo_free(matcoo* m);
matcoo* matcoo_zeroes(int dx, int dy);
matcoo* matcoo_build(matcoo* m, float val, int i, int j);
float matcoo_get(matcoo* m, int mi, int mj);
void matcoo_print(matcoo* m);

#endif
