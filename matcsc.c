#include "matcsc.h"


matcsc *matcsc_new(const float *data, int dimX, int dimY) {
    matcsc *m = (matcsc *) malloc(sizeof(matcsc));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dimX + dimY) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    list_append(&m->ia, 0.0f); // An extra element ia[0] = 0 is used by convention
    m->dimX = dimX;
    m->dimY = dimY;
    int stride = dimX;
    int ia = 0;
    for (int j = 0; j < dimX; j++) {
        for (int i = 0; i < dimY; i++) {
            float val = data[j * stride + i];
            if (val != 0) {
                list_append(&m->nnz, val);
                list_append(&m->ja, (float) i);
                ia++;
            }
        }
        list_append(&m->ia, (float) ia);
    }
    return m;
}

void matcsc_free(matcsc *m) {
    list_free(m->nnz);
    list_free(m->ia);
    list_free(m->ja);
    free(m);
}

matcsc *matcsc_zeroes(int dx, int dy) {
    matcsc *m = (matcsc *) malloc(sizeof(matcsc));
    int capacity = (dx + dy) / 2; // set initial capacity to average of col and row lengths
    int inc = (int) fmax((dx + dy) / 2.0, 5.0); // set increment to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    for (int i = 0; i < dx + 1; i++) { // An extra element ia[0] = 0 is used by convention, hence + 1
        list_append(&m->ia, 0.0f);
    }
    m->dimX = dx;
    m->dimY = dy;
    return m;
}

float matcsc_get(matcsc *m, int i, int j) {
    int ia = (int) list_get(m->ia, j); // the number of values in nnz to skip to get to col i
    int ia_next = (int) list_get(m->ia, j + 1); // the number of values in nnz to skip to get to col i + 1
    int els_this_col = ia_next - ia;
    if (els_this_col == 0) {
        return 0.0f;
    }
    for (int k = ia; k < ia_next; k++) {
        if ((int) list_get(m->ja, k) == i) {
            return list_get(m->nnz, k);
        }
    }
    return 0.0f;
}

void matcsc_rawprint(matcsc *m) {
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

void matcsc_print(matcsc *m) {
    for (int j = 0; j < m->dimX; j++) {
        for (int i = 0; i < m->dimY; i++) {
            printf("%.2f\t\t", matcsc_get(m, i, j));
        }
        printf("\n");
    }
}

matcsc *matcsc_sm(matcsc *m, float s) {
    for (int i = 0; i < m->nnz->length; i++) {
        int pi = i / m->nnz->increment; // part index, int division floors result
        int si = i - (m->nnz->increment * pi); // sub index, the index within the part
        m->nnz->parts[pi][si] *= s;
    }
    return m;
}

float matcsc_trace(matcsc *m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    float trace = 0;
    for (int i = 0; i < m->dimY; i++) {
        trace += matcsc_get(m, i, i);
    }
    return trace;
}

matcsc *matcsc_add(matcsc *m1, matcsc *m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {  // matrices must be of same rank to be added
        assert("Add is only defined for two matrices of the same rank.\n" && 0);
    }

    matcsc *mres = matcsc_zeroes(m1->dimX, m1->dimY);
    for (int j = 0; j < m1->dimX; j++) {
        // TODO: dont need to check rows where ia is same as last ia
        for (int i = 0; i < m1->dimY; i++) {
            float v1 = matcsc_get(m1, i, j);
            float v2 = matcsc_get(m1, i, j);
            float r = v1 + v2;
            if (r != 0) {
                list_append(&mres->nnz, r); // add the non zero value to nnz
                list_append(&mres->ja, (float) i); // add the row index to ja
            }
        }
        list_set(mres->ia, j + 1, (float) mres->nnz->length); // add the number of values to far to ia
    }
    return mres;
}

matcsc *matcsc_transpose(matcsc *m) {
    assert("NYI" && 0);
}

matcsc *matcsc_multiply(matcsc *m1, matcsc *m2) {
    assert("NYI" && 0);
}
