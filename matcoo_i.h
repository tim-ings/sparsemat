#ifndef SPARSEMAT_MATCOO_I_H
#define SPARSEMAT_MATCOO_I_H

#include "ll_int.h"
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>


typedef struct matcoo_i_ matcoo_i;
struct matcoo_i_ {
    ll_int *vals;
    ll_int *is;
    ll_int *js;
    int dimX;
    int dimY;
};

matcoo_i *matcoo_i_new(const int *data, int dimX, int dimY);

void matcoo_i_free(matcoo_i *m);

matcoo_i *matcoo_i_zeroes(int dx, int dy);

matcoo_i *matcoo_i_build(matcoo_i *m, int val, int i, int j);

int matcoo_i_get(matcoo_i *m, int mi, int mj);

void matcoo_i_print(matcoo_i *m);

matcoo_i *matcoo_i_sm(matcoo_i *m, float a);

int matcoo_i_trace(matcoo_i *m);

matcoo_i *matcoo_i_add(matcoo_i *m1, matcoo_i *m2);

matcoo_i *matcoo_i_transpose(matcoo_i *m);

matcoo_i *matcoo_i_multiply(matcoo_i *m1, matcoo_i *m2);

bool matcoo_i_equals(matcoo_i *m1, matcoo_i *m2);

#endif
