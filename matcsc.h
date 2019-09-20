#ifndef SPARSEMAT_MATCSC_H
#define SPARSEMAT_MATCSC_H
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include "ll_float.h"

typedef struct matcsc_ matcsc;

struct matcsc_ {
    float* nnz; // The non-zero values stored in column wise ordering order (left to right, top to bottom)
    float* ia; // The number of elements in each column. An extra element IA[0] = 0 is used by convention. This array can be used to index into the NNZ array for each i-th column
    float* ja; // Stores the row index of each non-zero element
    int nnz_len;
    int ia_len;
    int ja_len;
    int dimX;
    int dimY;
};

matcsc* matcsc_new(const float* data, int dimX, int dimY);
matcsc* matcsc_sm(matcsc* m, float s);
void matcsc_print(matcsc* m);
float matcsc_trace(matcsc* m);

#endif
