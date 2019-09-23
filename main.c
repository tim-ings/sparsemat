#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "matcoo.h"
#include "matcsr.h"
#include "matcsc.h"
#include <omp.h>


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
        for (int j = 0; j < mcoo->dimX; j++) {
            float v1 = matcoo_get(mcoo, i, j);
            float v2 = matcsr_get(mcsr, i, j, 0);
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
    for (int j = 0; j < mcoo->dimX; j++) {
        for (int i = 0; i < mcoo->dimY; i++) {
            float v1 = matcoo_get(mcoo, i, j);
            float v2 = matcsc_get(mcsc, j, i, 0);
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

matcsr *matcsr_fromfile(const char *fp) {
    int dx, dy;
    float *d = readfile(fp, &dx, &dy);
    matcsr *m = matcsr_new(d, dx, dy);
    return m;
}

matcsc *matcsc_fromfile(const char *fp) {
    int dx, dy;
    float *d = readfile(fp, &dx, &dy);
    matcsc *m = matcsc_new(d, dx, dy);
    return m;
}

matcoo *mat_multiply(matcsr *m1, matcsc *m2) {
    if (m1->dimX != m2->dimY) {
        assert("Multiplication is only defined for two matrices where A.dimX == B.dimY" && 0);
    }
    matcoo *mres = matcoo_zeroes(m2->dimX, m1->dimY);
    int i;
#pragma omp parallel for num_threads(16) private(i) shared(mres, m1, m2) default(none) collapse(2)
    for (i = 0; i < mres->dimY; i++) {
        for (int j = 0; j < mres->dimX; j++) {
            float sum = 0.0f;
            int found1 = 0;
            int found2 = 0;
            for (int k = 0; k < m1->dimX; k++) {
                float v1 = matcsr_get(m1, i, k, &found1);
                float v2 = matcsc_get(m2, j, k, &found2);
                sum += v1 * v2;
            }
            matcoo_build(mres, sum, i, j);
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

int test_all(const char *f1, const char *f2) {
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

int time_csr(const char *f1, const char *f2) {
    // load
    timer_start();
    matcsr *m_init = matcsr_fromfile(f1);
    printf("CSR Load: %.4f\n", timer_end());
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

    // transpose
//    matcsr *m_transpose = matcsr_fromfile(f1);
//    timer_start();
//    matcsr *m_transpose_res = matcsr_transpose(m_transpose);
//    printf("CSR transpose: %.4f\n", timer_end());
//    matcsr_free(m_transpose_res);

    // multiply
//    matcsr *m_multiply1 = matcsr_fromfile(f1);
//    matcsr *m_multiply2 = matcsr_fromfile(f2);
//    timer_start();
//    matcsr *m_multiply_res = matcsr_multiply(m_multiply1, m_multiply2);
//    printf("CSR multiply: %.4f\n", timer_end());
//    matcsr_free(m_multiply_res);

    matcsr *mm_csr = matcsr_fromfile(f1);
    matcsc *mm_csc = matcsc_fromfile(f2);
    timer_start();
    matcoo *mm_res = mat_multiply(mm_csr, mm_csc);
    printf("CSR x CSC multiply: %.4f\n", timer_end());
//    matcsr_free(mm_csr);
//    matcsc_free(mm_csc);
    return 0;
}

int time_csc(const char *f1, const char *f2) {
    // load
    timer_start();
    matcsc *m_init = matcsc_fromfile(f1);
    printf("CSC Load: %.4f\n", timer_end());
    matcsc_free(m_init);

    // scalar multiplication
    matcsc *m_sm = matcsc_fromfile(f1);
    timer_start();
    matcsc *m_sm_res = matcsc_sm(m_sm, 2);
    printf("CSC sm: %.4f\n", timer_end());
    matcsc_free(m_sm_res);

    // trace
    matcsc *m_trace = matcsc_fromfile(f1);
    timer_start();
    matcsc_trace(m_trace);
    printf("CSC trace: %.4f\n", timer_end());
    matcsc_free(m_trace);

    // add
    matcsc *m_add1 = matcsc_fromfile(f1);
    matcsc *m_add2 = matcsc_fromfile(f2);
    timer_start();
    matcsc *m_add_res = matcsc_add(m_add1, m_add2);
    printf("CSC add: %.4f\n", timer_end());
    matcsc_free(m_add_res);
    matcsc_free(m_add2);

    // transpose
//    matcsc *m_transpose = matcsc_fromfile(f1);
//    timer_start();
//    matcsc *m_transpose_res = matcsc_transpose(m_transpose);
//    printf("CSC transpose: %.4f\n", timer_end());
//    matcsc_free(m_transpose_res);

    // multiply
//    matcsc *m_multiply1 = matcsc_fromfile(f1);
//    matcsc *m_multiply2 = matcsc_fromfile(f2);
//    timer_start();
//    matcsc *m_multiply_res = matcsc_multiply(m_multiply1, m_multiply2);
//    printf("CSC multiply: %.4f\n", timer_end());
//    matcsc_free(m_multiply_res);
    return 0;
}

int main(int argc, char *argv[]) {
//    const char *f1 = "data/float64.in";
//    const char *f2 = "data/float64.in";
//    test_all(f1, f2);

//    time_coo("data/float64.in", "data/float64.in");
//    printf("\n");
    time_csr("data/float256.in", "data/float256.in");
//    printf("\n");
//    time_csc("data/float128.in", "data/float128.in");

    return EXIT_SUCCESS;
}