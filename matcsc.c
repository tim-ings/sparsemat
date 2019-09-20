#include "matcsc.h"


matcsc* matcsc_new(const float* data, int dimX, int dimY) {
    matcsc* m = (matcsc*)malloc(sizeof(matcsc));
    m->dimX = dimX;
    m->dimY = dimY;
    m->nnz = malloc(sizeof(float) * dimX * dimY); // allocate a worst case array
    m->ia = malloc(sizeof(float) * dimX * dimY); // allocate a worst case array
    m->ja = malloc(sizeof(float) * dimX * dimY); // allocate a worst case array
    m->nnz_len = 0;
    m->ia_len = 0;
    m->ja_len = 0;
    int stride = dimX;
    for (int i = 0; i < dimX; i++) {
        for (int j = 0; j < dimY; j++) {
            float val = data[i * stride + j];
            if (val != 0) {
                m->nnz[m->nnz_len++] = val;
                m->ia[m->ia_len++] = (float)j;
            }
        }
        m->ja[m->ja_len++] = (float)m->nnz_len;
    }
    return m;
}

void matcsc_print(matcsc* m) {
    int next_val_i = 0;
    int els_so_far_i = 1;
    int next_row_index_i = 0;
    float next_val = m->nnz[next_val_i++];
    int els_so_far = (int)m->ia[els_so_far_i++];
    int next_row_index = (int)m->ja[next_row_index_i++];
    int vals_this_col = 0;
    for (int j = 0; j < m->dimX; j++) {
        for (int i = 0; i < m->dimY; i++) {
            if (vals_this_col < els_so_far && i == next_row_index) {
                printf("%.2f\t\t", next_val);
                next_val = m->nnz[next_val_i++];
                next_row_index = (int)m->ja[next_row_index_i++];
                vals_this_col++;
            } else {
                printf("0.00\t\t");
            }
        }
        els_so_far = (int)m->ia[els_so_far_i++];
        printf("\n");
    }
}

matcsc* matcsc_sm(matcsc* m, float s) {
    for (int i = 0; i < m->nnz_len; i++) {
        m->nnz[i] *= s;
    }
    return m;
}

float matcsc_trace(matcsc* m) {
    float sum = 0;
    int next_val_i = 0;
    int els_so_far_i = 1;
    int next_row_index_i = 0;
    float next_val = m->nnz[next_val_i++];
    int els_so_far = (int)m->ia[els_so_far_i++];
    int next_row_index = (int)m->ja[next_row_index_i++];
    int vals_this_col = 0;
    for (int j = 0; j < m->dimX; j++) {
        for (int i = 0; i < m->dimY; i++) {
            if (vals_this_col < els_so_far && i == next_row_index) {
                if (i == j) {
                    sum += next_val;
                }
                next_val = m->nnz[next_val_i++];
                next_row_index = (int)m->ja[next_row_index_i++];
                vals_this_col++;
            }
        }
        els_so_far = (int)m->ia[els_so_far_i++];
    }
    return sum;
}
