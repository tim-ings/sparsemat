#ifndef SPARSEMAT_MATCSC_I_H
#define SPARSEMAT_MATCSC_I_H

#include <stdlib.h>
#include <stdio.h>
#include "ll_int.h"
#include "list_i.h"
#include <math.h>
#include <assert.h>

typedef struct matcsc_i_ matcsc_i;

struct matcsc_i_ {
    list_i *nnz; // The non-zero values stored in column wise ordering order (top to bottom, left to right)
    list_i *ia; // number of elements in nnz, used to index into nnz array for each i-th col
    list_i *ja; // nnz row indices
    int dimX;
    int dimY;
};

matcsc_i *matcsc_i_new(const int *data, int dimX, int dimY);

void matcsc_i_free(matcsc_i *m);

matcsc_i *matcsc_i_zeroes(int dx, int dy);

int matcsc_i_get(matcsc_i *m, int i, int j, int* col_skip);

void matcsc_i_rawprint(matcsc_i *m);

void matcsc_i_print(matcsc_i *m);

matcsc_i *matcsc_i_sm(matcsc_i *m, float s);

int matcsc_i_trace(matcsc_i *m);

matcsc_i *matcsc_i_add(matcsc_i *m1, matcsc_i *m2);

matcsc_i *matcsc_i_transpose(matcsc_i *m);

matcsc_i *matcsc_i_multiply(matcsc_i *m1, matcsc_i *m2);

#endif
