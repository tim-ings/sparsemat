#include "matcsc_i.h"


matcsc_i *matcsc_i_new(const int *data, int dimX, int dimY) {
    matcsc_i *m = (matcsc_i *) malloc(sizeof(matcsc_i));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_i_new(capacity, inc);
    m->ia = list_i_new(capacity, inc);
    m->ja = list_i_new(capacity, inc);
    list_i_append(&m->ia, 0); // An extra element ia[0] = 0 is used by convention
    m->dimX = dimX;
    m->dimY = dimY;
    int stride = dimX;
    int ia = 0;
    for (int j = 0; j < dimX; j++) {
        for (int i = 0; i < dimY; i++) {
            int val = data[j * stride + i];
            if (val != 0) {
                list_i_append(&m->nnz, val);
                list_i_append(&m->ja, i);
                ia++;
            }
        }
        list_i_append(&m->ia, ia);
    }
    return m;
}

void matcsc_i_free(matcsc_i *m) {
    list_i_free(m->nnz);
    list_i_free(m->ia);
    list_i_free(m->ja);
    free(m);
}

matcsc_i *matcsc_i_zeroes(int dx, int dy) {
    matcsc_i *m = (matcsc_i *) malloc(sizeof(matcsc_i));
    int capacity = (dx + dy) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dx + dy) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_i_new(capacity, inc);
    m->ia = list_i_new(capacity, inc);
    m->ja = list_i_new(capacity, inc);
    for (int i = 0; i < dx + 1; i++) { // An extra element ia[0] = 0 is used by convention, hence + 1
        list_i_append(&m->ia, 0);
    }
    m->dimX = dx;
    m->dimY = dy;
    return m;
}

int matcsc_i_get(matcsc_i *m, int i, int j, int *col_skip) {
    int ia = list_i_get(m->ia, j); // the number of values in nnz to skip to get to col i
    int ia_next = list_i_get(m->ia, j + 1); // the number of values in nnz to skip to get to col i + 1
    if (ia_next - ia == 0)
        return 0; // there are no nnz's this row
    if (col_skip)
        ia += *col_skip;
    for (int k = ia; k < ia_next; k++) {
        if (list_i_get(m->ja, k) == i) {
            if (col_skip)
                *col_skip += 1;
            return list_i_get(m->nnz, k);
        }
    }
    return 0;
}

void matcsc_i_rawprint(matcsc_i *m) {
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

void matcsc_i_print(matcsc_i *m) {
    for (int j = 0; j < m->dimX; j++) {
        for (int i = 0; i < m->dimY; i++) {
            printf("%.2f\t\t", matcsc_i_get(m, i, j, 0));
        }
        printf("\n");
    }
}

matcsc_i *matcsc_i_sm(matcsc_i *m, float s) {
    for (int i = 0; i < m->nnz->length; i++) {
        int pi = i / m->nnz->increment; // part index, int division floors result
        int si = i - (m->nnz->increment * pi); // sub index, the index within the part
        m->nnz->parts[pi][si] *= s;
    }
    return m;
}

int matcsc_i_trace(matcsc_i *m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    int trace = 0;
    for (int i = 0; i < m->dimY; i++) {
        trace += matcsc_i_get(m, i, i, 0);
    }
    return trace;
}

matcsc_i *matcsc_i_add(matcsc_i *m1, matcsc_i *m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {  // matrices must be of same rank to be added
        assert("Add is only defined for two matrices of the same rank.\n" && 0);
    }

    matcsc_i *mres = matcsc_i_zeroes(m1->dimX, m1->dimY);
    for (int j = 0; j < m1->dimX; j++) {
        // TODO: dont need to check rows where ia is same as last ia
        int found1, found2;
        for (int i = 0; i < m1->dimY; i++) {
            int v1 = matcsc_i_get(m1, i, j, &found1);
            int v2 = matcsc_i_get(m1, i, j, &found2);
            int r = v1 + v2;
            if (r != 0) {
                list_i_append(&mres->nnz, r); // add the non zero value to nnz
                list_i_append(&mres->ja, i); // add the row index to ja
            }
        }
        list_i_set(mres->ia, j + 1, mres->nnz->length); // add the number of values to far to ia
    }
    return mres;
}

matcsc_i *matcsc_i_transpose(matcsc_i *m) {
    assert("NYI" && 0);
}

matcsc_i *matcsc_i_multiply(matcsc_i *m1, matcsc_i *m2) {
    assert("NYI" && 0);
}
