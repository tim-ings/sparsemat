#include "matcsr.h"


matcsr *matcsr_new(const float *data, int dimX, int dimY) {
    matcsr *m = (matcsr *) malloc(sizeof(matcsr));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    list_append(&m->ia, 0.0f); // An extra element ia[0] = 0 is used by convention
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
        list_append(&m->ia, (float) ia);
    }
    return m;
}

void matcsr_free(matcsr *m) {
    list_free(m->nnz);
    list_free(m->ia);
    list_free(m->ja);
    free(m);
}

matcsr *matcsr_zeroes(int dimX, int dimY) {
    matcsr *m = (matcsr *) malloc(sizeof(matcsr));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    for (int i = 0; i < dimY + 1; i++) { // An extra element ia[0] = 0 is used by convention, hence + 1
        list_append(&m->ia, 0.0f);
    }
    m->dimX = dimX;
    m->dimY = dimY;
    return m;
}

float matcsr_get(matcsr *m, int i, int j, int *row_skip) {
    int ia = (int) list_get(m->ia, i); // the number of values in nnz to skip to get to row i
    int ia_next = (int) list_get(m->ia, i + 1); // the number of values in nnz to skip to get to row i + 1
    if (ia_next - ia == 0)
        return 0.0f; // there are no nnz's this row
    if (row_skip)
        ia += *row_skip;
    for (int k = ia; k < ia_next; k++) {
        if ((int) list_get(m->ja, k) == j) {
            if (row_skip)
                *row_skip += 1;
            return list_get(m->nnz, k);
        }
    }
    return 0.0f;
}

void matcsr_rawprint(matcsr *m) {
    printf("nnz \t= ");
    for (int i = 0; i < m->nnz->length; i++) {
        float v = list_get(m->nnz, i);
        printf("%.2f \t", v);
    }
    printf("\n");
    printf("ia[%d]\t= ", m->ia->length);
    for (int i = 0; i < m->ia->length; i++) {
        float v = list_get(m->ia, i);
        printf("%d    \t", (int) v);
    }
    printf("\n");
    printf("ja[%d]\t= ", m->ja->length);
    for (int i = 0; i < m->ja->length; i++) {
        float v = list_get(m->ja, i);
        printf("%d    \t", (int) v);
    }
    printf("\n");
}

void matcsr_print(matcsr *m) {
    for (int i = 0; i < m->dimY; i++) {
        int found = 0;
        for (int j = 0; j < m->dimX; j++) {
            printf("%.2f\t\t", matcsr_get(m, i, j, &found));
        }
        printf("\n");
    }
}

matcsr *matcsr_sm(matcsr *m, float s) {
    for (int i = 0; i < m->nnz->length; i++) {
        int pi = i / m->nnz->increment; // part index, int division floors result
        int si = i - (m->nnz->increment * pi); // sub index, the index within the part
        m->nnz->parts[pi][si] *= s;
    }
    return m;
}

float matcsr_trace(matcsr *m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    float trace = 0;
    int f;
    for (int i = 0; i < m->dimY; i++) {
        trace += matcsr_get(m, i, i, &f);
    }
    return trace;
}

matcsr *matcsr_add(matcsr *m1, matcsr *m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {  // matrices must be of same rank to be added
        assert("Add is only defined for two matrices of the same rank.\n" && 0);
    }

    matcsr *mres = matcsr_zeroes(m1->dimX, m1->dimY);
    for (int i = 0; i < m1->dimY; i++) {
        // TODO: dont need to check rows where ia is same as last ia
        int found1, found2;
        for (int j = 0; j < m1->dimX; j++) {
            float v1 = matcsr_get(m1, i, j, &found1);
            float v2 = matcsr_get(m1, i, j, &found2);
            float r = v1 + v2;
            if (r != 0) {
                list_append(&mres->nnz, r); // add the non zero value to nnz
                list_append(&mres->ja, (float) j); // add the col index to ja
            }
        }
        list_set(mres->ia, i + 1, (float) mres->nnz->length); // add the number of values to far to ia
    }
    return mres;
}

matcsr *matcsr_transpose(matcsr *m) {
    assert("NYI" && 0);
}

matcsr *matcsr_multiply(matcsr *m1, matcsr *m2) {
    assert("NYI" && 0);
}