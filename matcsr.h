#ifndef SPARSEMAT_MATCSR_H
#define SPARSEMAT_MATCSR_H
#include <stdlib.h>
#include <stdio.h>
#include "ll_float.h"

typedef struct matcsr_ matcsr;

struct matcsr_ {
    float* nnz; // The non-zero values stored in row-major
                // order (left to right, top to bottom)
    float* ia;  // The number of elements in each row. An extra element IA[0] = 0 is
                // used by convention. This array can be used to index into the NNZ array for each i-th row
    float* ja;  // Stores the column index of each non-zero element
    int dimX;
    int dimY;
};

matcsr* matcsr_new(const float* data, int dimX, int dimY);
matcsr* matcsr_sm(matcsr* m, float s);
void matcsr_print(matcsr* m);


#endif //SPARSEMAT_MATCSR_H
