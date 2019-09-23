#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "matcoo.h"
#include "matcoo_i.h"
#include "matcsr.h"
#include "matcsr_i.h"
#include "matcsc.h"
#include "matcsc_i.h"


struct timespec tstart = {0, 0}, tend = {0, 0};

void timer_start() {
    clock_gettime(CLOCK_MONOTONIC, &tstart);
}

double timer_end() {
    clock_gettime(CLOCK_MONOTONIC, &tend);
    return ((double) tend.tv_sec + 1.0e-9 * (double) tend.tv_nsec) -
           ((double) tstart.tv_sec + 1.0e-9 * (double) tstart.tv_nsec);
}

bool compare_coo_csr(matcoo *mcoo, matcsr *mcsr) {
    if (mcoo->dimX != mcsr->dimX || mcoo->dimY != mcsr->dimY) {
        return false; // mats must have same rank to be equal
    }
    // check all the elements
    for (int i = 0; i < mcoo->dimY; i++) {
        int f2 = 0;
        for (int j = 0; j < mcoo->dimX; j++) {
            float v1 = matcoo_get(mcoo, i, j);
            float v2 = matcsr_get(mcsr, i, j, &f2);
            if (v1 != v2) {
                return false;
            }
        }
    }
    return true;
}

bool compare_coo_csr_i(matcoo_i *mcoo, matcsr_i *mcsr) {
    if (mcoo->dimX != mcsr->dimX || mcoo->dimY != mcsr->dimY) {
        return false; // mats must have same rank to be equal
    }
    // check all the elements
    for (int i = 0; i < mcoo->dimY; i++) {
        int f2 = 0;
        for (int j = 0; j < mcoo->dimX; j++) {
            int v1 = matcoo_i_get(mcoo, i, j);
            int v2 = matcsr_i_get(mcsr, i, j, &f2);
            if (v1 != v2) {
                return false;
            }
        }
    }
    return true;
}

bool compare_coo_csc(matcoo *mcoo, matcsc *mcsc) {
    if (mcoo->dimX != mcsc->dimX || mcoo->dimY != mcsc->dimY) {
        return false; // mats must have same rank to be equal
    }
    // check all the elements
    for (int i = 0; i < mcoo->dimY; i++) {
        int f2 = 0;
        for (int j = 0; j < mcoo->dimX; j++) {
            float v1 = matcoo_get(mcoo, i, j);
            float v2 = matcsc_get(mcsc, j, i, &f2);
            if (v1 != v2) {
                return false;
            }
        }
    }
    return true;
}

bool compare_coo_csc_i(matcoo_i *mcoo, matcsc_i *mcsc) {
    if (mcoo->dimX != mcsc->dimX || mcoo->dimY != mcsc->dimY) {
        return false; // mats must have same rank to be equal
    }
    // check all the elements
    for (int j = 0; j < mcoo->dimX; j++) {
        int f2 = 0;
        for (int i = 0; i < mcoo->dimY; i++) {
            int v1 = matcoo_i_get(mcoo, i, j);
            int v2 = matcsc_i_get(mcsc, j, i, &f2);
            if (v1 != v2) {
                return false;
            }
        }
    }
    return true;
}

float *readfile(const char *fileName, int *out_dimX, int *out_dimY) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s\n", fileName);
        return NULL;
    }
    int counter = 0;
    float *data = NULL;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[read - 1] = '\0';
        switch (counter) {
            case 0: {
                break;
            }
            case 1: {
                *out_dimX = (int) strtol(line, NULL, 10);
                break;
            }
            case 2: {
                *out_dimY = (int) strtol(line, NULL, 10);
                data = (float *) malloc(sizeof(float) * (*out_dimX) * (*out_dimY));
                break;
            }
            case 3: {
                const char delim[2] = " ";
                char *token = strtok(line, delim);
                char *err;
                data[0] = (float) strtod(token, &err);
                for (int i = 1; i < ((*out_dimX) * (*out_dimY)); i++) {
                    token = strtok(NULL, delim);
                    data[i] = (float) strtod(token, &err);
                }
                break;
            }
            default: {
                printf("Unexpected line in %s, ignoring it\n", fileName);
                return NULL;
            }
        }
        counter++;
    }
    fclose(fp);
    if (line)
        free(line);

    return data;
}

int *readfile_i(const char *fileName, int *out_dimX, int *out_dimY) {
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE *fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s\n", fileName);
        return NULL;
    }
    int counter = 0;
    int *data = NULL;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[read - 1] = '\0';
        switch (counter) {
            case 0: {
                break;
            }
            case 1: {
                *out_dimX = (int) strtol(line, NULL, 10);
                break;
            }
            case 2: {
                *out_dimY = (int) strtol(line, NULL, 10);
                data = malloc(sizeof(int) * (*out_dimX) * (*out_dimY));
                break;
            }
            case 3: {
                const char delim[2] = " ";
                char *token = strtok(line, delim);
                char *err;
                data[0] = (float) strtod(token, &err);
                for (int i = 1; i < ((*out_dimX) * (*out_dimY)); i++) {
                    token = strtok(NULL, delim);
                    data[i] = (float) strtol(token, &err, 10);
                }
                break;
            }
            default: {
                printf("Unexpected line in %s, ignoring it\n", fileName);
                return NULL;
            }
        }
        counter++;
    }
    fclose(fp);
    if (line)
        free(line);

    return data;
}

matcoo *matcoo_fromfile(const char *fp) {
    int dx, dy;
    float *d = readfile(fp, &dx, &dy);
    matcoo *m = matcoo_zeroes(dx, dy);

    for (int i = 0; i < dy; i++) {
        for (int j = 0; j < dx; j++) {
            float v = d[i * dx + j];
            m = matcoo_build(m, v, i, j);
        }
    }
    return m;
}

matcoo_i *matcoo_i_fromfile(const char *fp) {
    int dx, dy;
    int *d = readfile_i(fp, &dx, &dy);
    matcoo_i *m = matcoo_i_zeroes(dx, dy);

    for (int i = 0; i < dy; i++) {
        for (int j = 0; j < dx; j++) {
            int v = d[i * dx + j];
            m = matcoo_i_build(m, v, i, j);
        }
    }
    return m;
}

matcsr *matcsr_fromfile(const char *fp) {
    int dx, dy;
    float *d = readfile(fp, &dx, &dy);
    matcsr *m = matcsr_new(d, dx, dy);
    return m;
}

matcsr_i *matcsr_i_fromfile(const char *fp) {
    int dx, dy;
    int *d = readfile_i(fp, &dx, &dy);
    matcsr_i *m = matcsr_i_new(d, dx, dy);
    return m;
}

matcsc *matcsc_fromfile(const char *fp) {
    int dx, dy;
    float *d = readfile(fp, &dx, &dy);
    matcsc *m = matcsc_new(d, dx, dy);
    return m;
}

matcsc_i *matcsc_i_fromfile(const char *fp) {
    int dx, dy;
    int *d = readfile_i(fp, &dx, &dy);
    matcsc_i *m = matcsc_i_new(d, dx, dy);
    return m;
}

matcoo *mat_multiply(matcsr *m1, matcsc *m2) {
    if (m1->dimX != m2->dimY) {
        assert("Multiplication is only defined for two matrices where A.dimX == B.dimY" && 0);
    }
    matcoo *mres = matcoo_zeroes(m2->dimX, m1->dimY);
    int i;
//#pragma omp parallel for num_threads(16) private(i, found1, found2) shared(mres, m1, m2) default(none) collapse(2)
    for (i = 0; i < mres->dimY; i++) {
        for (int j = 0; j < mres->dimX; j++) {
            float sum = 0.0f;
            for (int k = 0; k < m1->dimX; k++) {
                float v1 = matcsr_get(m1, i, k, 0);
                float v2 = matcsc_get(m2, j, k, 0);
                sum += v1 * v2;
            }
            matcoo_build(mres, sum, i, j);
        }
    }
    return mres;
}

matcoo_i *mat_multiply_i(matcsr_i *m1, matcsc_i *m2) {
    if (m1->dimX != m2->dimY) {
        assert("Multiplication is only defined for two matrices where A.dimX == B.dimY" && 0);
    }
    matcoo_i *mres = matcoo_i_zeroes(m2->dimX, m1->dimY);
    int i;
#pragma omp parallel for num_threads(16) private(i) shared(mres, m1, m2) default(none) collapse(2)
    for (i = 0; i < mres->dimY; i++) {
        for (int j = 0; j < mres->dimX; j++) {
            int sum = 0;
            int found1 = 0;
            int found2 = 0;
            for (int k = 0; k < m1->dimX; k++) {
                int v1 = matcsr_i_get(m1, i, k, &found1);
                int v2 = matcsc_i_get(m2, j, k, &found2);
                sum += v1 * v2;
            }
            matcoo_i_build(mres, sum, i, j);
        }
    }
    return mres;
}

int test_coo(const char *f1, const char *f2) {
    printf("\tTest COO: Load:\n");
    matcoo *m_init = matcoo_fromfile(f1);
    matcoo_print(m_init);
    matcoo_free(m_init);

    printf("\tTest COO: Scalar Multiplication:\n");
    matcoo *m_sm = matcoo_fromfile(f1);
    matcoo *m_sm_res = matcoo_sm(m_sm, 2);
    matcoo_print(m_sm_res);
    matcoo_free(m_sm_res);

    printf("\tTest COO: Trace:\n");
    matcoo *m_trace = matcoo_fromfile(f1);
    float t = matcoo_trace(m_trace);
    printf("t=%.2f\n", t);
    matcoo_free(m_trace);

    printf("\tTest COO: Add:\n");
    matcoo *m_add1 = matcoo_fromfile(f1);
    matcoo *m_add2 = matcoo_fromfile(f2);
    matcoo *m_add_res = matcoo_add(m_add1, m_add2);
    matcoo_print(m_add_res);
    matcoo_free(m_add_res);
    matcoo_free(m_add2);

    printf("\tTest COO: Transpose:\n");
    matcoo *m_transpose = matcoo_fromfile(f1);
    matcoo *m_transpose_res = matcoo_transpose(m_transpose);
    matcoo_print(m_transpose_res);
    matcoo_free(m_transpose_res);

    printf("\tTest COO: Multiply:\n");
    matcoo *m_multiply1 = matcoo_fromfile(f1);
    matcoo *m_multiply2 = matcoo_fromfile(f2);
    matcoo *m_multiply_res = matcoo_multiply(m_multiply1, m_multiply2);
    matcoo_print(m_multiply_res);
    matcoo_free(m_multiply_res);
    return 0;
}

int test_csr(const char *f1, const char *f2) {
    printf("\tTest CSR: Load:\n");
    matcsr *m_init = matcsr_fromfile(f1);
    matcsr_print(m_init);
    matcsr_free(m_init);

    printf("\tTest CSR: Scalar Multiplication:\n");
    matcsr *m_sm = matcsr_fromfile(f1);
    matcsr *m_sm_res = matcsr_sm(m_sm, 2);
    matcsr_print(m_sm_res);
    matcsr_free(m_sm_res);

    printf("\tTest CSR: Trace:\n");
    matcsr *m_trace = matcsr_fromfile(f1);
    float t = matcsr_trace(m_trace);
    printf("t=%.2f\n", t);
    matcsr_free(m_trace);

    printf("\tTest CSR: Add:\n");
    matcsr *m_add1 = matcsr_fromfile(f1);
    matcsr *m_add2 = matcsr_fromfile(f2);
    matcsr *m_add_res = matcsr_add(m_add1, m_add2);
    matcsr_print(m_add_res);
    matcsr_free(m_add_res);
    matcsr_free(m_add2);

//    printf("\tTest CSR: Transpose:\n");
//    matcsr *m_transpose = matcsr_fromfile(f1);
//    matcsr *m_transpose_res = matcsr_transpose(m_transpose);
//    matcsr_print(m_transpose_res);
//    matcsr_free(m_transpose_res);

//    printf("\tTest CSR: Multiply:\n");
//    matcsr *m_multiply1 = matcsr_fromfile(f1);
//    matcsr *m_multiply2 = matcsr_fromfile(f2);
//    matcsr *m_multiply_res = matcsr_multiply(m_multiply1, m_multiply2);
//    matcsr_print(m_multiply_res);
//    matcsr_free(m_multiply_res);
    return 0;
}

int test_csc(const char *f1, const char *f2) {
    printf("\tTest CSC: Load:\n");
    matcsc *m_init = matcsc_fromfile(f1);
    matcsc_print(m_init);
    matcsc_free(m_init);

    printf("\tTest CSC: Scalar Multiplication:\n");
    matcsc *m_sm = matcsc_fromfile(f1);
    matcsc *m_sm_res = matcsc_sm(m_sm, 2);
    matcsc_print(m_sm_res);
    matcsc_free(m_sm_res);

    printf("\tTest CSC: Trace:\n");
    matcsc *m_trace = matcsc_fromfile(f1);
    float t = matcsc_trace(m_trace);
    printf("t=%.2f\n", t);
    matcsc_free(m_trace);

    printf("\tTest CSC: Add:\n");
    matcsc *m_add1 = matcsc_fromfile(f1);
    matcsc *m_add2 = matcsc_fromfile(f2);
    matcsc *m_add_res = matcsc_add(m_add1, m_add2);
    matcsc_print(m_add_res);
    matcsc_free(m_add_res);
    matcsc_free(m_add2);

//    printf("\tTest CSC: Transpose:\n");
//    matcsc *m_transpose = matcsc_fromfile(f1);
//    matcsc *m_transpose_res = matcsc_transpose(m_transpose);
//    matcsc_print(m_transpose_res);
//    matcsc_free(m_transpose_res);
//
//    printf("\tTest CSC: Multiply:\n");
//    matcsc *m_multiply1 = matcsc_fromfile(f1);
//    matcsc *m_multiply2 = matcsc_fromfile(f2);
//    matcsc *m_multiply_res = matcsc_multiply(m_multiply1, m_multiply2);
//    matcsc_print(m_multiply_res);
//    matcsc_free(m_multiply_res);
    return 0;
}

int test_all_float(const char *f1, const char *f2) {
    // loading
    matcoo *init_coo = matcoo_fromfile(f1);
    matcsr *init_csr = matcsr_fromfile(f1);
    matcsc *init_csc = matcsc_fromfile(f1);
    assert(compare_coo_csr(init_coo, init_csr));
    assert(compare_coo_csc(init_coo, init_csc));
    matcoo_free(init_coo);
    matcsr_free(init_csr);
    matcsc_free(init_csc);

    // scalar multiplication
    float a = 2.0f;
    matcoo *sm_coo = matcoo_fromfile(f1);
    matcsr *sm_csr = matcsr_fromfile(f1);
    matcsc *sm_csc = matcsc_fromfile(f1);
    sm_coo = matcoo_sm(sm_coo, a);
    sm_csr = matcsr_sm(sm_csr, a);
    sm_csc = matcsc_sm(sm_csc, a);
    assert(compare_coo_csr(sm_coo, sm_csr));
    assert(compare_coo_csc(sm_coo, sm_csc));
    matcoo_free(sm_coo);
    matcsr_free(sm_csr);
    matcsc_free(sm_csc);

    // trace
    matcoo *trace_coo = matcoo_fromfile(f1);
    matcsr *trace_csr = matcsr_fromfile(f1);
    matcsc *trace_csc = matcsc_fromfile(f1);
    float coo_trace = matcoo_trace(trace_coo);
    float csr_trace = matcsr_trace(trace_csr);
    float csc_trace = matcsc_trace(trace_csc);
    assert(coo_trace == csr_trace);
    assert(coo_trace == csc_trace);
    matcoo_free(trace_coo);
    matcsr_free(trace_csr);
    matcsc_free(trace_csc);

    // add
    matcoo *coo_add1 = matcoo_fromfile(f1);
    matcoo *coo_add2 = matcoo_fromfile(f2);
    matcoo *coo_add_res = matcoo_add(coo_add1, coo_add2);
    matcsr *csr_add1 = matcsr_fromfile(f1);
    matcsr *csr_add2 = matcsr_fromfile(f2);
    matcsr *csr_add_res = matcsr_add(csr_add1, csr_add2);
    matcsc *csc_add1 = matcsc_fromfile(f1);
    matcsc *csc_add2 = matcsc_fromfile(f2);
    matcsc *csc_add_res = matcsc_add(csc_add1, csc_add2);
    assert(compare_coo_csr(coo_add_res, csr_add_res));
    assert(compare_coo_csc(coo_add_res, csc_add_res));
    matcoo_free(coo_add1);
    matcoo_free(coo_add2);
    matcsr_free(csr_add1);
    matcsr_free(csr_add2);
    matcsc_free(csc_add1);
    matcsc_free(csc_add2);

    // transpose
    matcoo *coo_transpose = matcoo_fromfile(f1);
    matcoo *coo_transpose_res = matcoo_transpose(coo_transpose);
//    matcsr *csr_transpose = matcsr_fromfile(f1);
//    matcsr *csr_transpose_res = matcsr_transpose(csr_transpose);
//    assert(compare_coo_csr(coo_transpose_res, csr_transpose_res));
    matcoo_free(coo_transpose);
//    matcsr_free(csr_transpose);

    // multiply
    matcoo *coo_multi1 = matcoo_fromfile(f1);
    matcoo *coo_multi2 = matcoo_fromfile(f2);
    matcoo *coo_multi_res = matcoo_multiply(coo_multi1, coo_multi2);
    matcsr *csr_multi1 = matcsr_fromfile(f1);
    matcsc *csc_multi2 = matcsc_fromfile(f2);
    matcoo *coo_multi_res2 = mat_multiply(csr_multi1, csc_multi2);
    assert(matcoo_equals(coo_multi_res, coo_multi_res2));
    matcoo_free(coo_multi1);
    matcoo_free(coo_multi2);
    matcsr_free(csr_multi1);
    matcsc_free(csc_multi2);
    matcoo_free(coo_multi_res);
    matcoo_free(coo_multi_res2);
    return 0;
}

int test_all_int(const char *f1, const char *f2) {
    // loading
    matcoo_i *init_coo = matcoo_i_fromfile(f1);
    matcsr_i *init_csr = matcsr_i_fromfile(f1);
    matcsc_i *init_csc = matcsc_i_fromfile(f1);
    assert(compare_coo_csr_i(init_coo, init_csr));
    assert(compare_coo_csc_i(init_coo, init_csc));
    matcoo_i_free(init_coo);
    matcsr_i_free(init_csr);
    matcsc_i_free(init_csc);

    // scalar multiplication
    float a = 2.0f;
    matcoo_i *sm_coo = matcoo_i_fromfile(f1);
    matcsr_i *sm_csr = matcsr_i_fromfile(f1);
    matcsc_i *sm_csc = matcsc_i_fromfile(f1);
    sm_coo = matcoo_i_sm(sm_coo, a);
    sm_csr = matcsr_i_sm(sm_csr, a);
    sm_csc = matcsc_i_sm(sm_csc, a);
    assert(compare_coo_csr_i(sm_coo, sm_csr));
    assert(compare_coo_csc_i(sm_coo, sm_csc));
    matcoo_i_free(sm_coo);
    matcsr_i_free(sm_csr);
    matcsc_i_free(sm_csc);

    // trace
    matcoo_i *trace_coo = matcoo_i_fromfile(f1);
    matcsr_i *trace_csr = matcsr_i_fromfile(f1);
    matcsc_i *trace_csc = matcsc_i_fromfile(f1);
    int coo_trace = matcoo_i_trace(trace_coo);
    int csr_trace = matcsr_i_trace(trace_csr);
    int csc_trace = matcsc_i_trace(trace_csc);
    assert(coo_trace == csr_trace);
    assert(coo_trace == csc_trace);
    matcoo_i_free(trace_coo);
    matcsr_i_free(trace_csr);
    matcsc_i_free(trace_csc);

    // add
    matcoo_i *coo_add1 = matcoo_i_fromfile(f1);
    matcoo_i *coo_add2 = matcoo_i_fromfile(f2);
    matcoo_i *coo_add_res = matcoo_i_add(coo_add1, coo_add2);
    matcsr_i *csr_add1 = matcsr_i_fromfile(f1);
    matcsr_i *csr_add2 = matcsr_i_fromfile(f2);
    matcsr_i *csr_add_res = matcsr_i_add(csr_add1, csr_add2);
    matcsc_i *csc_add1 = matcsc_i_fromfile(f1);
    matcsc_i *csc_add2 = matcsc_i_fromfile(f2);
    matcsc_i *csc_add_res = matcsc_i_add(csc_add1, csc_add2);
    assert(compare_coo_csr_i(coo_add_res, csr_add_res));
    assert(compare_coo_csc_i(coo_add_res, csc_add_res));
    matcoo_i_free(coo_add1);
    matcoo_i_free(coo_add2);
    matcsr_i_free(csr_add1);
    matcsr_i_free(csr_add2);
    matcsc_i_free(csc_add1);
    matcsc_i_free(csc_add2);

    // transpose
    matcoo_i *coo_transpose = matcoo_i_fromfile(f1);
    matcoo_i *coo_transpose_res = matcoo_i_transpose(coo_transpose);
//    matcsr *csr_transpose = matcsr_fromfile(f1);
//    matcsr *csr_transpose_res = matcsr_transpose(csr_transpose);
//    assert(compare_coo_csr(coo_transpose_res, csr_transpose_res));
    matcoo_i_free(coo_transpose);
//    matcsr_free(csr_transpose);

    // multiply
    matcoo_i *coo_multi1 = matcoo_i_fromfile(f1);
    matcoo_i *coo_multi2 = matcoo_i_fromfile(f2);
    matcoo_i *coo_multi_res = matcoo_i_multiply(coo_multi1, coo_multi2);
    matcsr_i *csr_multi1 = matcsr_i_fromfile(f1);
    matcsc_i *csc_multi2 = matcsc_i_fromfile(f2);
    matcoo_i *coo_multi_res2 = mat_multiply_i(csr_multi1, csc_multi2);
    assert(matcoo_i_equals(coo_multi_res, coo_multi_res2));
    matcoo_i_free(coo_multi1);
    matcoo_i_free(coo_multi2);
    matcsr_i_free(csr_multi1);
    matcsc_i_free(csc_multi2);
    matcoo_i_free(coo_multi_res);
    matcoo_i_free(coo_multi_res2);
    return 0;
}

int time_coo(const char *f1, const char *f2) {
    // load
    timer_start();
    matcoo *m_init = matcoo_fromfile(f1);
    printf("COO Load: %.4f\n", timer_end());
    matcoo_free(m_init);

    // scalar multiplication
    matcoo *m_sm = matcoo_fromfile(f1);
    timer_start();
    matcoo *m_sm_res = matcoo_sm(m_sm, 2);
    printf("COO sm: %.4f\n", timer_end());
    matcoo_free(m_sm_res);

    // trace
    matcoo *m_trace = matcoo_fromfile(f1);
    timer_start();
    matcoo_trace(m_trace);
    printf("COO trace: %.4f\n", timer_end());
    matcoo_free(m_trace);

    // add
    matcoo *m_add1 = matcoo_fromfile(f1);
    matcoo *m_add2 = matcoo_fromfile(f2);
    timer_start();
    matcoo *m_add_res = matcoo_add(m_add1, m_add2);
    printf("COO add: %.4f\n", timer_end());
    matcoo_free(m_add_res);
    matcoo_free(m_add2);

    // transpose
    matcoo *m_transpose = matcoo_fromfile(f1);
    timer_start();
    matcoo *m_transpose_res = matcoo_transpose(m_transpose);
    printf("COO transpose: %.4f\n", timer_end());
    matcoo_free(m_transpose_res);

    // multiply
    matcoo *m_multiply1 = matcoo_fromfile(f1);
    matcoo *m_multiply2 = matcoo_fromfile(f2);
    timer_start();
    matcoo *m_multiply_res = matcoo_multiply(m_multiply1, m_multiply2);
    printf("COO multiply: %.4f\n", timer_end());
    matcoo_free(m_multiply_res);
    return 0;
}

int time_csr_float(const char *f1, const char *f2) {
    // load
    timer_start();
    matcsr *m_init = matcsr_fromfile(f1);
    printf("CSR Load: %.4f\n", timer_end());
//    matcsr_print(m_init);
//    matcsr_free(m_init);

    // scalar multiplication
    matcsr *m_sm = matcsr_fromfile(f1);
    timer_start();
    matcsr *m_sm_res = matcsr_sm(m_sm, 2);
    printf("CSR sm: %.4f\n", timer_end());
//    matcsr_free(m_sm_res);

    // trace
    matcsr *m_trace = matcsr_fromfile(f1);
    timer_start();
    matcsr_trace(m_trace);
    printf("CSR trace: %.4f\n", timer_end());
//    matcsr_free(m_trace);

    // add
    matcsr *m_add1 = matcsr_fromfile(f1);
    matcsr *m_add2 = matcsr_fromfile(f2);
    timer_start();
    matcsr *m_add_res = matcsr_add(m_add1, m_add2);
    printf("CSR add: %.4f\n", timer_end());
//    matcsr_free(m_add_res);
//    matcsr_free(m_add2);

    matcsr *mm_csr = matcsr_fromfile(f1);
    matcsc *mm_csc = matcsc_fromfile(f2);
    timer_start();
    matcoo *mm_res = mat_multiply(mm_csr, mm_csc);
    printf("CSR x CSC multiply: %.4f\n", timer_end());
//    matcsr_free(mm_csr);
//    matcsc_free(mm_csc);


    printf("CSR FLOAT PASSED ALL TESTS\n");
    return 0;
}

int time_csr_int(const char *f1, const char *f2) {
    // load
    timer_start();
    matcsr_i *m_init = matcsr_i_fromfile(f1);
    printf("CSR INT Load: %.4f\n", timer_end());
//    matcsr_i_print(m_init);
//    matcsr_free(m_init);

    // scalar multiplication
    matcsr_i *m_sm = matcsr_i_fromfile(f1);
    timer_start();
    matcsr_i *m_sm_res = matcsr_i_sm(m_sm, 2);
    printf("CSR INT sm: %.4f\n", timer_end());
//    matcsr_i_free(m_sm_res);

    // trace
    matcsr_i *m_trace = matcsr_i_fromfile(f1);
    timer_start();
    matcsr_i_trace(m_trace);
    printf("CSR INT trace: %.4f\n", timer_end());
//    matcsr_i_free(m_trace);

    // add
    matcsr_i *m_add1 = matcsr_i_fromfile(f1);
    matcsr_i *m_add2 = matcsr_i_fromfile(f2);
    timer_start();
    matcsr_i *m_add_res = matcsr_i_add(m_add1, m_add2);
    printf("CSR INT add: %.4f\n", timer_end());
//    matcsr_i_free(m_add_res);
//    matcsr_i_free(m_add2);

    matcsr_i *mm_csr = matcsr_i_fromfile(f1);
    matcsc_i *mm_csc = matcsc_i_fromfile(f2);
    timer_start();
    matcoo_i *mm_res = mat_multiply_i(mm_csr, mm_csc);
    printf("CSR INT x CSC multiply: %.4f\n", timer_end());
//    matcsr_i_free(mm_csr);
//    matcsc_i_free(mm_csc);

    printf("CSR INT PASSED ALL TESTS\n");
    return 0;
}

bool is_float(char* filepath) {
    char line[10];
    FILE *fp = fopen(filepath, "r");
    if (fp) {
        fgets(line, 10, fp);
        printf("Data from the file:\n%s", line);
        fclose(fp);
    }
}

int run_sm(char* filepath) {
    if (is_float(filepath)) {
        matcoo *m = matcoo_fromfile(filepath);
        timer_start();
        m = matcoo_sm(m, 2);
        printf("Operation completed in %.4fs\n", timer_end());
        matcoo_print(m);
        matcoo_free(m);
    } else {
        matcoo_i *m = matcoo_i_fromfile(filepath);
        timer_start();
        m = matcoo_i_sm(m, 2);
        printf("Operation completed in %.4fs\n", timer_end());
        matcoo_i_print(m);
//        matcsr_i_free(m);
    }
}

int print_usage() {
    printf("Usage:\n"
           "\t--sm \tScalar multiplication\n"
           "\t--tr \tTrace\n"
           "\t--ad \tAddition\n"
           "\t--ts \tTranspose\n"
           "\t--mm \tMatrix multiplication\n");
}

int main(int argc, char *argv[]) {
    char operation = '0';
    int file_count = 1;
    char* file1 = NULL;
    char* file2 = NULL;
    // scalar multiplication
    if (argc < 2) {
        print_usage();
    } else {
        // go over every arg
        for (int i = 0; i < argc; i++) {
            char* v = argv[i];
            if (strcmp(v, "--sm") == 0) {
                operation = 's';
            } else if (strcmp(v, "--tr") == 0) {
                operation = 'r';
            } else if (strcmp(v, "--ad") == 0) {
                operation = 'a';
                file_count = 2;
            } else if (strcmp(v, "--ts") == 0) {
                operation = 't';
            } else if (strcmp(v, "--mm") == 0) {
                operation = 'm';
                file_count = 2;
            } else if (strcmp(v, "-f") == 0) {
                file1 = argv[i + 1];
                if (file_count > 1) {
                    file2 = argv[i + 2];
                }
            }
        }
    }
    if (file_count == 2) {
        printf("Running operation %c on %d files: %s and %s\n", operation, file_count, file1, file2);
    } else {
        printf("Running operation %c on %d files: %s\n", operation, file_count, file1);
    }
    switch (operation) {
        case 's':
            run_sm(file1);
            break;
        default:
            break;
    }


//    const char *f1 = "data/float123.in";
//    const char *f2 = "data/float123.in";
//    test_all_float(f1, f2);
//
//    time_csr_float("data/float256.in", "data/float256.in");
//    time_csr_int("data/int256.in", "data/int256.in");

    return EXIT_SUCCESS;
}