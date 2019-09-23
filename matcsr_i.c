#include "matcsr_i.h"


matcsr_i *matcsr_i_new(const int *data, int dimX, int dimY) {
    matcsr_i *m = (matcsr_i *) malloc(sizeof(matcsr_i));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_i_new(capacity, inc);
    m->ia = list_i_new(capacity, inc);
    m->ja = list_i_new(capacity, inc);
    list_i_append(&m->ia, 0); // An extra element ia[0] = 0 is used by convention
    m->dimX = dimX;
    m->dimY = dimY;
    int stride = dimY;
    int ia = 0;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            int val = data[i * stride + j];
            if (val != 0) {
                list_i_append(&m->nnz, val);
                list_i_append(&m->ja, j);
                ia++;
            }
        }
        list_i_append(&m->ia, ia);
    }
    return m;
}

void matcsr_i_free(matcsr_i *m) {
    list_i_free(m->nnz);
    list_i_free(m->ia);
    list_i_free(m->ja);
    free(m);
}

matcsr_i *matcsr_i_zeroes(int dimX, int dimY) {
    matcsr_i *m = (matcsr_i *) malloc(sizeof(matcsr_i));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_i_new(capacity, inc);
    m->ia = list_i_new(capacity, inc);
    m->ja = list_i_new(capacity, inc);
    for (int i = 0; i < dimY + 1; i++) { // An extra element ia[0] = 0 is used by convention, hence + 1
        list_i_append(&m->ia, 0);
    }
    m->dimX = dimX;
    m->dimY = dimY;
    return m;
}

int matcsr_i_get(matcsr_i *m, int i, int j, int *row_skip) {
    int ia = list_i_get(m->ia, i); // the number of values in nnz to skip to get to row i
    int ia_next = list_i_get(m->ia, i + 1); // the number of values in nnz to skip to get to row i + 1
    if (ia_next - ia == 0)
        return 0; // there are no nnz's this row
    if (row_skip)
        ia += *row_skip;
    for (int k = ia; k < ia_next; k++) {
        if (list_i_get(m->ja, k) == j) {
            if (row_skip)
                *row_skip += 1;
            return list_i_get(m->nnz, k);
        }
    }
    return 0;
}

void matcsr_i_rawprint(matcsr_i *m) {
    printf("nnz \t= ");
    for (int i = 0; i < m->nnz->length; i++) {
        int v = list_i_get(m->nnz, i);
        printf("%.2f \t", v);
    }
    printf("\n");
    printf("ia[%d]\t= ", m->ia->length);
    for (int i = 0; i < m->ia->length; i++) {
        int v = list_i_get(m->ia, i);
        printf("%d    \t", (int) v);
    }
    printf("\n");
    printf("ja[%d]\t= ", m->ja->length);
    for (int i = 0; i < m->ja->length; i++) {
        int v = list_i_get(m->ja, i);
        printf("%d    \t", (int) v);
    }
    printf("\n");
}

void matcsr_i_print(matcsr_i *m) {
    for (int i = 0; i < m->dimY; i++) {
        int found = 0;
        for (int j = 0; j < m->dimX; j++) {
            printf("%d\t\t", matcsr_i_get(m, i, j, &found));
        }
        printf("\n");
    }
}

matcsr_i *matcsr_i_sm(matcsr_i *m, float s) {
    for (int i = 0; i < m->nnz->length; i++) {
        int pi = i / m->nnz->increment; // part index, int division floors result
        int si = i - (m->nnz->increment * pi); // sub index, the index within the part
        m->nnz->parts[pi][si] *= s;
    }
    return m;
}

int matcsr_i_trace(matcsr_i *m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    int trace = 0;
    int f;
    for (int i = 0; i < m->dimY; i++) {
        trace += matcsr_i_get(m, i, i, &f);
    }
    return trace;
}

matcsr_i *matcsr_i_add(matcsr_i *m1, matcsr_i *m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {  // matrices must be of same rank to be added
        assert("Add is only defined for two matrices of the same rank.\n" && 0);
    }

    matcsr_i *mres = matcsr_i_zeroes(m1->dimX, m1->dimY);
    for (int i = 0; i < m1->dimY; i++) {
        // TODO: dont need to check rows where ia is same as last ia
        int found1, found2;
        for (int j = 0; j < m1->dimX; j++) {
            int v1 = matcsr_i_get(m1, i, j, &found1);
            int v2 = matcsr_i_get(m1, i, j, &found2);
            int r = v1 + v2;
            if (r != 0) {
                list_i_append(&mres->nnz, r); // add the non zero value to nnz
                list_i_append(&mres->ja, j); // add the col index to ja
            }
        }
        list_i_set(mres->ia, i + 1, mres->nnz->length); // add the number of values to far to ia
    }
    return mres;
}

matcsr_i *matcsr_i_transpose(matcsr_i *m) {
    assert("NYI" && 0);
}

matcsr_i *matcsr_i_multiply(matcsr_i *m1, matcsr_i *m2) {
    assert("NYI" && 0);
}