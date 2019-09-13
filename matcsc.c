#include "matcsc.h"


matcsc* matcsc_new(const float* data, int dimX, int dimY) {
#ifdef DEBUG_VERBOSE
    printf("\nBuilding new %dx%d CSC matrix\n", dimX, dimY);
#endif
    // build a linked list of all the non zero values (for nnz)
    ll_float* nzlist = ll_float_new();
    // build a list of the number of elements per column so far (for ia)
    ll_float* elements_per_col = ll_float_new();
    // build a list of the row index of each non zero value (for ja)
    ll_float* col_indices = ll_float_new();
    int stride = dimY;
    int els_so_far = 0;
    int prints_so_far = 0;
    for (int j = 0; j < dimX; j++) {
        for (int i = 0; i < dimY; i++) {
            float val = data[i * stride + j];
            if (val != 0) {
                ll_float_push(nzlist, val);
                ll_float_push(col_indices, (float)j);
                els_so_far++;
                if (prints_so_far < PRINT_COUNT_MAX / 4) {
#ifdef DEBUG_VERBOSE
                    printf("{ nzv: %.2f, els_so_far: %d, col_index: %d }, ", val, els_so_far, j);
#endif
                    prints_so_far++;
                }
            }
        }
        ll_float_push(elements_per_col, (float)els_so_far);
    }
#ifdef DEBUG_VERBOSE
    printf("...\n");
#endif

    matcsc* m = (matcsc*)malloc(sizeof(matcsc));
    m->dimX = dimX;
    m->dimY = dimY;

    // transform linked list into array
    m->nnz = malloc(sizeof(float) * nzlist->length);
    int i = 0;
    prints_so_far = 0;
#ifdef DEBUG_VERBOSE
    printf("nnz: [ ");
#endif
    ll_float_node* nnz_cur = nzlist->first;
    while (nnz_cur) {
        float val = ll_float_next(&nnz_cur);
        m->nnz[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
#ifdef DEBUG_VERBOSE
            printf("%.2f, ", val);
#endif
            prints_so_far++;
        }
    }
#ifdef DEBUG_VERBOSE
    printf("]...\n");
#endif

    // transform linked list into array
    m->ia = malloc(sizeof(float) * elements_per_col->length + 1);
    m->ia[0] = 0;
    i = 1;
    prints_so_far = 0;
#ifdef DEBUG_VERBOSE
    printf("ia:  [ ");
    printf("%.2f, ", m->ia[0]);
#endif
    ll_float_node* epr_cur = elements_per_col->first;
    while (epr_cur) {
        float val = ll_float_next(&epr_cur);
        m->ia[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
#ifdef DEBUG_VERBOSE
            printf("%.2f, ", val);
#endif
            prints_so_far++;
        }
    }
#ifdef DEBUG_VERBOSE
    printf("]...\n");
#endif

    // transform linked list into array
    m->ja = malloc(sizeof(float) * col_indices->length);
    i = 0;
    prints_so_far = 0;
#ifdef DEBUG_VERBOSE
    printf("ja:  [ ");
#endif
    ll_float_node* ci_cur = col_indices->first;
    while (ci_cur) {
        float val = ll_float_next(&ci_cur);
        m->ja[i] = val;
        i++;
        if (prints_so_far < PRINT_COUNT_MAX) {
#ifdef DEBUG_VERBOSE
            printf("%.2f, ", val);
#endif
            prints_so_far++;
        }
    }

#ifdef DEBUG_VERBOSE
    printf("]...\n%dx%d CSC matrix built\n", dimX, dimY);
#endif
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
