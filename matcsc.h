#ifndef SPARSEMAT_MATCSC_H
#define SPARSEMAT_MATCSC_H

#include <stdlib.h>
#include <stdio.h>
#include "ll_float.h"
#include "list.h"
#include <math.h>
#include <assert.h>

typedef struct matcsc_ matcsc;

struct matcsc_ {
    list *nnz; // The non-zero values stored in column wise ordering order (top to bottom, left to right)
    list *ia; // number of elements in nnz, used to index into nnz array for each i-th col
    list *ja; // nnz row indices
    int dimX;
    int dimY;
};

matcsc *matcsc_new(const float *data, int dimX, int dimY);

void matcsc_free(matcsc *m);

matcsc *matcsc_zeroes(int dx, int dy);

float matcsc_get(matcsc *m, int i, int j, int* col_skip);

void matcsc_rawprint(matcsc *m);

void matcsc_print(matcsc *m);

matcsc *matcsc_sm(matcsc *m, float s);

float matcsc_trace(matcsc *m);

matcsc *matcsc_add(matcsc *m1, matcsc *m2);

matcsc *matcsc_transpose(matcsc *m);

matcsc *matcsc_multiply(matcsc *m1, matcsc *m2);

#endif
