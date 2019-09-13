#include "matcsr.h"


matcsr* matcsr_new(const float* data, int dimX, int dimY) {
    printf("\nBuilding new %dx%d CSR matrix\n", dimX, dimY);
    // build a linked list of all the non zero values (for nnz)
    ll_float* nzlist = ll_float_new();
    // build a list of the number of elements per row so far (for ia)
    ll_float* elements_per_row = ll_float_new();
    // build a list of the column index of each non zero value (for ja)
    ll_float* col_indices = ll_float_new();
    int stride = dimY;
    int els_so_far = 0;
    int prints_so_far = 0;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float val = data[i * stride + j];
            if (val != 0) {
                ll_float_push(nzlist, val);
                ll_float_push(col_indices, (float)j);
                els_so_far++;
                if (prints_so_far < PRINT_COUNT_MAX / 4) {
                    printf("{ nzv: %f, els_so_far: %d, col_index: %d }, ", val, els_so_far, j);
                    prints_so_far++;
                }
            }
        }
        ll_float_push(elements_per_row, (float)els_so_far);
    }
    printf("...\n");

    matcsr* m = (matcsr*)malloc(sizeof(matcsr));

    // transform linked list into array
    m->nnz = malloc(sizeof(float) * nzlist->length);
    int i = 0;
    prints_so_far = 0;
    printf("nnz: [ ");
    ll_float_node* nnz_cur = nzlist->first;
    while (nnz_cur) {
        float val = ll_float_next(&nnz_cur);
        m->nnz[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
            printf("%.2f, ", val);
            prints_so_far++;
        }
    }
    printf("]...\n");

    // transform linked list into array
    m->ia = malloc(sizeof(float) * elements_per_row->length + 1);
    m->ia[0] = 0;
    i = 1;
    prints_so_far = 0;
    printf("ia:  [ ");
    printf("%.2f, ", m->ia[0]);
    ll_float_node* epr_cur = elements_per_row->first;
    while (epr_cur) {
        float val = ll_float_next(&epr_cur);
        m->ia[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
            printf("%.2f, ", val);
            prints_so_far++;
        }
    }
    printf("]...\n");

    // transform linked list into array
    m->ja = malloc(sizeof(float) * col_indices->length);
    i = 0;
    prints_so_far = 0;
    printf("ja:  [ ");
    ll_float_node* ci_cur = col_indices->first;
    while (ci_cur) {
        float val = ll_float_next(&ci_cur);
        m->ja[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
            printf("%.2f, ", val);
            prints_so_far++;
        }
    }

    printf("]...\n%dx%d CSR matrix built\n", dimX, dimY);
    return m;
}
