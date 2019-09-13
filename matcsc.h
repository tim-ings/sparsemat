#ifndef SPARSEMAT_MATCSC_H
#define SPARSEMAT_MATCSC_H
#include <stdlib.h>
#include <stdio.h>
#include "ll_float.h"

typedef struct matcsc_ matcsc;

struct matcsc_ {
    float* nnz; // The non-zero values stored in column wise ordering
                // order (left to right, top to bottom)
    float* ia;  // The number of elements in each column. An extra element IA[0] = 0 is
                // used by convention. This array can be used to index into the NNZ array for each i-th column
    float* ja;  // Stores the row index of each non-zero element
    int dimX;
    int dimY;
};

matcsc* matcsc_new(const float* data, int dimX, int dimY);
void matcsc_print(matcsc* m);

#endif //SPARSEMAT_MATCSC_H
