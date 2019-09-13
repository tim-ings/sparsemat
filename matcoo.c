#include "matcoo.h"


matcoo* matcoo_new(const float* data, int dimX, int dimY) {
    printf("\nBuilding new %dx%d COO matrix\n[ ", dimX, dimY);
    ll_float* coords_val = ll_float_new();
    ll_float* coords_i = ll_float_new();
    ll_float* coords_j = ll_float_new();
    int stride = dimY;
    int prints_so_far = 0;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float val = data[i * stride + j];
            if (val != 0) {
                ll_float_push(coords_val, val);
                ll_float_push(coords_i, (float)i);
                ll_float_push(coords_j, (float)j);
                if (prints_so_far < PRINT_COUNT_MAX / 1.5) {
                    printf("(%d, %d, %.2f), ", j, i, val);
                    prints_so_far++;
                }
            }
        }
    }
    ll_float_push(coords_val, data[dimY * stride + dimX]);
    ll_float_push(coords_i, (float)dimY);
    ll_float_push(coords_j, (float)dimX);

    matcoo* m = (matcoo*)malloc(sizeof(matcoo));
    m->dimX = dimX;
    m->dimY = dimY;
    m->coords_val = coords_val;
    m->coords_i = coords_i;
    m->coords_j = coords_j;

    printf("]...\n%dx%d COO matrix built\n", dimX, dimY);
    return m;
}
