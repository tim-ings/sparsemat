#include "matcsr.h"


matcsr* matcsr_new(const float* data, int dimX, int dimY) {
    matcsr* m = (matcsr*)malloc(sizeof(matcsr));
    int size = (dimX + dimY) / 2;
    int inc = (int)fmax((dimX + dimY) / 4.0, 5.0);
    m->nnz = list_new(size, inc);
    m->ia = list_new(size, inc); // this call returns garbage
    m->ja = list_new(size, inc);
    m->ia = list_new(size, inc); // so we have to call it again after ja
    list_append(&m->ia, 0.0f); // An extra element ia[0] = 0 is used by convention
    list* myl = list_new(size, inc);
    m->dimX = dimX;
    m->dimY = dimY;
    int stride = dimY;
    int ia = 0;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float val = data[i * stride + j];
            if (val != 0) {
                list_append(&m->nnz, val);
                list_append(&m->ja, (float) j);
                ia++;
            }
        }
        list_append(&myl, (float) ia);
    }
    m->ia = myl;
    for (int i = 0; i < m->nnz->capacity; i++) {
        printf("nnz[%d] = %.2f\n", i, list_get(m->nnz, i));
    }
    return m;
}

void matcsr_print(matcsr* m) {
    int sofar = 0;
    for (int i = 0; i < m->dimY; i++) {
        int row_offset = (int) list_get(m->ia, i); // the number of values in nnz to skip to get to row i
        for (int j = 0; j < m->dimX; j++) {
            int col_ind = (int) list_get(m->ja, sofar); // column index of the next non 0 value
            if (sofar < row_offset && j == col_ind) { // is this i,j a non zero value?
                float v = list_get(m->nnz, sofar++); // get the value and increment sofar
                printf("%.2f\t\t", v);
            } else { // the value is not non-0
                printf("0.00\t\t");
            }
        }
        printf("\n");
    }

//    int next_val_i = 0;
//    int els_so_far_i = 1;
//    int next_col_index_i = 0;
//    float next_val = m->nnz[next_val_i++];
//    int els_so_far = (int)m->ia[els_so_far_i++];
//    int next_col_index = (int)m->ja[next_col_index_i++];
//    int vals_this_row = 0;
//    for (int i = 0; i < m->dimY; i++) {
//        for (int j = 0; j < m->dimX; j++) {
//            if (vals_this_row < els_so_far && j == next_col_index) {
//                printf("%.2f\t\t", next_val);
//                next_val = m->nnz[next_val_i++];
//                next_col_index = (int)m->ja[next_col_index_i++];
//                vals_this_row++;
//            } else {
//                printf("0.00\t\t");
//            }
//        }
//        els_so_far = (int)m->ia[els_so_far_i++];
//        printf("\n");
//    }
}