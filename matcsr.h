#ifndef SPARSEMAT_MATCSR_H
#define SPARSEMAT_MATCSR_H

#include <stdlib.h>
#include <stdio.h>
#include "ll_float.h"
#include "list.h"
#include <math.h>
#include <assert.h>


typedef struct matcsr_ matcsr;

struct matcsr_ {
    list *nnz; // The non-zero values stored in row-major order (left to right, top to bottom)
    list *ia;   // The number of elements in each row. An extra element IA[0] = 0 is used by convention.
    // This array can be used to index into the NNZ array for each i-th row
    list *ja; // Stores the column index of each non-zero element
    int dimX;
    int dimY;
};

matcsr *matcsr_new(const float *data, int dimX, int dimY);

void matcsr_free(matcsr *m);

matcsr *matcsr_zeroes(int dx, int dy);

float matcsr_get(matcsr *m, int mi, int mj);

void matcsr_rawprint(matcsr *m);

void matcsr_print(matcsr *m);

matcsr *matcsr_sm(matcsr *m, float s);

float matcsr_trace(matcsr *m);

matcsr *matcsr_add(matcsr *m1, matcsr *m2);

matcsr *matcsr_transpose(matcsr *m);

matcsr *matcsr_multiply(matcsr *m1, matcsr *m2);

#endif
