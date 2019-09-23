#ifndef SPARSEMAT_MATCSR_I_H
#define SPARSEMAT_MATCSR_I_H

#include <stdlib.h>
#include <stdio.h>
#include "list_i.h"
#include <math.h>
#include <assert.h>
#include <stdbool.h>


typedef struct matcsr_i_ matcsr_i;

struct matcsr_i_ {
    list_i *nnz; // The non-zero values stored in row-major order (left to right, top to bottom)
    list_i *ia;   // The number of elements in each row. An extra element IA[0] = 0 is used by convention.
    // This array can be used to index into the NNZ array for each i-th row
    list_i *ja; // Stores the column index of each non-zero element
    int dimX;
    int dimY;
};

matcsr_i *matcsr_i_new(const int *data, int dimX, int dimY);

void matcsr_i_free(matcsr_i *m);

matcsr_i *matcsr_i_zeroes(int dx, int dy);

int matcsr_i_get(matcsr_i *m, int i, int j, int *row_skip);

void matcsr_i_rawprint(matcsr_i *m);

void matcsr_i_print(matcsr_i *m);

matcsr_i *matcsr_i_sm(matcsr_i *m, float s);

int matcsr_i_trace(matcsr_i *m);

matcsr_i *matcsr_i_add(matcsr_i *m1, matcsr_i *m2);

matcsr_i *matcsr_i_transpose(matcsr_i *m);

matcsr_i *matcsr_i_multiply(matcsr_i *m1, matcsr_i *m2);

#endif
