#include "matcsr.h"


matcsr* matcsr_fromfile(const char* filepath) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(filepath, "r");
    if (fp == NULL) {
        printf("File does not exist: %s", filepath);
        return NULL;
    }
    matcsr* m = NULL;
    int dimX = -1;
    int dimY = -1;
    int counter = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[read - 1] = '\0';
        switch (counter) {
            case 0: {
                break;
            }
            case 1: {
                dimX = (int)strtol(line, NULL, 10);
                break;
            }
            case 2: {
                dimY = (int)strtol(line, NULL, 10);
                m = matcsr_zeroes(dimX, dimY);
                break;
            }
            case 3: {
                const char delim[2] = " ";
                char* err;
                for (int i = 0; i < dimY; i++) {
                    for (int j = 0; j < dimX; j++) {
                        char* token = strtok(i == 0 && j == 0 ? line : NULL, delim);
                        float val = (float)strtod(token, &err);
                    }
                }
                break;
            }
            default: {
                printf("Unexpected line in %s, ignoring it\n", filepath);
                return NULL;
            }
        }
        counter++;
    }
    fclose(fp);
    if (line)
        free(line);
    return m;
}

matcsr* matcsr_new(const float* data, int dimX, int dimY) {
    matcsr* m = (matcsr*)malloc(sizeof(matcsr));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int)fmax((dimX + dimY) / 2.0, 5.0); // set initial capacity to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    list_set(m->ia, 0, 0.0f); // An extra element ia[0] = 0 is used by convention
    m->dimX = dimX;
    m->dimY = dimY;
    m->_ia_index = 0;
    int stride = dimY;
    int ia = 0;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            if (data) {
                float val = data[i * stride + j];
                if (val != 0) {
                    list_append(&m->nnz, val);
                    list_append(&m->ja, (float) j);
                    ia++;
                }
            }
        }
        list_append(&m->ia, (float) ia);
    }
    return m;
}

matcsr* matcsr_zeroes(int dimX, int dimY) {
    matcsr* m = (matcsr*)malloc(sizeof(matcsr));
    int capacity = (dimX + dimY) / 2; // set initial capacity to average of col and row lengths
    int inc = (int)fmax((dimX + dimY) / 2.0, 5.0); // set initial capacity to average of col and row lengths or 5
    m->nnz = list_new(capacity, inc);
    m->ia = list_new(capacity, inc);
    m->ja = list_new(capacity, inc);
    for (int i = 0; i < dimY + 1; i++) { // An extra element ia[0] = 0 is used by convention, hence + 1
        list_append(&m->ia, 0.0f);
    }
    m->dimX = dimX;
    m->dimY = dimY;
    m->_ia_index = 0;
    return m;
}

matcsr* matcsr_build(matcsr* m, float val, int i, int j) {
    if (val != 0) {
        list_set(m->ia, 0, 0.0f); // An extra element ia[0] = 0 is used by convention
        for (int ii = 0; ii < m->dimY; ii++) {
            for (int jj = 0; jj < m->dimX; jj++) {
                list_append(&m->nnz, val);
                list_append(&m->ja, (float) j);
                m->_ia_index++;
            }
        }
        list_append(&m->ia, (float) m->_ia_index);
    }
    return m;
}

void matcsr_rawprint(matcsr* m) {
    printf("nnz = ");
    for (int i = 0; i < m->nnz->length; i++) {
        float v = list_get(m->nnz, i);
        printf("%.2f \t", v);
    }
    printf("\n");
    printf(" ia[%d] = ", m->ia->length);
    for (int i = 0; i < m->ia->length; i++) {
        float v = list_get(m->ia, i);
        printf("%d    \t", (int)v);
    }
    printf("\n");
    printf(" ja[%d] = ", m->ja->length);
    for (int i = 0; i < m->ja->length; i++) {
        float v = list_get(m->ja, i);
        printf("%d    \t", (int)v);
    }
    printf("\n");

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
            } else { // the value is not non-0, so its 0
                printf("0.00\t\t");
            }
        }
        printf("\n");
    }
}

matcsr* matcsr_sm(matcsr* m, float s) {
    for (int i = 0; i < m->nnz->length; i++) {
        int pi = i / m->nnz->increment; // part index, int division floors result
        int si = i - (m->nnz->increment * pi); // sub index, the index within the part
        m->nnz->parts[pi][si] *= s;
    }
    return NULL;
}

matcsr* matcsr_add(matcsr* m1, matcsr* m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {
        return NULL; // matrices must be of same rank to be added
    }

    matcsr* mres = matcsr_new(NULL, m1->dimX, m1->dimY);

    int m1_sofar = 0;
    int m2_sofar = 0;
    for (int i = 0; i < m1->dimY; i++) {
        int m1_ro = (int) list_get(m1->ia, i); // the number of values in nnz to skip to get to row i in m1
        int m2_ro = (int) list_get(m2->ia, i); // the number of values in nnz to skip to get to row i in m2
        for (int j = 0; j < m1->dimX; j++) {
            int m1_col_ind = (int) list_get(m1->ja, m1_sofar); // column index of the next non 0 value
            int m2_col_ind = (int) list_get(m1->ja, m2_sofar); // column index of the next non 0 value
            float m1_val = 0;
            float m2_val = 0;
            // get the value from m1 at (i,j)
            if (m1_sofar < m1_ro && j == m1_col_ind) { // is this i,j a non zero value?
                m1_val = list_get(m1->nnz, m1_sofar++); // get the value and increment sofar
            }
            // get the value from m2 (i,j)
            if (m2_sofar < m2_ro && j == m2_col_ind) { // is this i,j a non zero value?
                m2_val = list_get(m2->nnz, m2_sofar++); // get the value and increment sofar
            }
//            mres = matscr_build(mres, i, j, m1_val + m2_val);
        }
        printf("\n");
    }
    return mres;
}

float matcsr_trace(matcsr* m) {
    if (m->dimX != m->dimY) {
        return -1; // matrices must be square for trace to be well defined
    }
    float trace = 0;
    int sofar = 0;
    for (int i = 0; i < m->dimY; i++) {
        int row_offset = (int) list_get(m->ia, i); // the number of values in nnz to skip to get to row i
        int col_ind = (int) list_get(m->ja, sofar); // column index of the next non 0 value
        while (col_ind < i) {
            col_ind = (int) list_get(m->ja, sofar++);
        }
        if (sofar < row_offset && i == col_ind) { // is this (i,j) a non zero value?
            float v = list_get(m->nnz, sofar++); // get the value and increment sofar
            trace += v;
        }
    }
    return trace;
}
