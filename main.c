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


#define STUDENT_NUMBER 21716194

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

matcoo *mat_multiply(matcsr *m1, matcsc *m2, int thread_count) {
    if (m1->dimX != m2->dimY) {
        assert("Multiplication is only defined for two matrices where A.dimX == B.dimY" && 0);
    }
    matcoo *mres = matcoo_zeroes(m2->dimX, m1->dimY);
    int i;
//#pragma omp parallel for num_threads(thread_count) private(i) shared(mres, m1, m2) default(none) collapse(2)
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

matcoo_i *mat_multiply_i(matcsr_i *m1, matcsc_i *m2, int thread_count) {
    if (m1->dimX != m2->dimY) {
        assert("Multiplication is only defined for two matrices where A.dimX == B.dimY" && 0);
    }
    matcoo_i *mres = matcoo_i_zeroes(m2->dimX, m1->dimY);
    int i;
    int threads = thread_count == -1 ? 8 : thread_count;
#pragma omp parallel for num_threads(threads) private(i) shared(mres, m1, m2) default(none) collapse(2)
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

// tests correctness of each matrix implementation as float
int test_correctness(const char *f1, const char *f2) {
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
    sm_csr = matcsr_sm(sm_csr, a, 8);
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
    matcoo_free(coo_transpose);

    // multiply
    matcoo *coo_multi1 = matcoo_fromfile(f1);
    matcoo *coo_multi2 = matcoo_fromfile(f2);
    matcoo *coo_multi_res = matcoo_multiply(coo_multi1, coo_multi2);
    matcsr *csr_multi1 = matcsr_fromfile(f1);
    matcsc *csc_multi2 = matcsc_fromfile(f2);
    matcoo *coo_multi_res2 = mat_multiply(csr_multi1, csc_multi2, -1);
    assert(matcoo_equals(coo_multi_res, coo_multi_res2));
    matcoo_free(coo_multi1);
    matcoo_free(coo_multi2);
    matcsr_free(csr_multi1);
    matcsc_free(csc_multi2);
    matcoo_free(coo_multi_res);
    matcoo_free(coo_multi_res2);
    return 0;
}

// tests correctness of each matrix implementation as int
int test_correctness_i(const char *f1, const char *f2) {
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
    matcoo_i_free(coo_transpose);

    // multiply
    matcoo_i *coo_multi1 = matcoo_i_fromfile(f1);
    matcoo_i *coo_multi2 = matcoo_i_fromfile(f2);
    matcoo_i *coo_multi_res = matcoo_i_multiply(coo_multi1, coo_multi2);
    matcsr_i *csr_multi1 = matcsr_i_fromfile(f1);
    matcsc_i *csc_multi2 = matcsc_i_fromfile(f2);
    matcoo_i *coo_multi_res2 = mat_multiply_i(csr_multi1, csc_multi2, -8);
    assert(matcoo_i_equals(coo_multi_res, coo_multi_res2));
    matcoo_i_free(coo_multi1);
    matcoo_i_free(coo_multi2);
    matcsr_i_free(csr_multi1);
    matcsc_i_free(csc_multi2);
    matcoo_i_free(coo_multi_res);
    matcoo_i_free(coo_multi_res2);
    return 0;
}

// checks if the given .in file is int or float
bool is_float(char *filepath) {
    char line[10];
    FILE *fp = fopen(filepath, "r");
    if (fp) {
        fgets(line, 10, fp);
        fclose(fp);
        return line[0] == 'f';
    }
}

char *get_log_name(const char* operation) {
    char *filename = malloc(sizeof(char) * 256);
    time_t t = time(NULL);
    struct tm mytm = *localtime(&t);
    sprintf(filename, "%d_%d%d%d_%d_%s.out", STUDENT_NUMBER, mytm.tm_mday, mytm.tm_mon + 1, mytm.tm_year + 1900,
            mytm.tm_sec * 1000, operation);
    return filename;
}

// runs the scalar multiplication operation
void run_sm(char *filepath, float a, bool should_print, bool should_log, int thread_count) {
    if (is_float(filepath)) {
        matcoo *m = matcoo_fromfile(filepath);
        timer_start();
        m = matcoo_sm(m, a);
        double elapsed = timer_end();
        printf("Scalar multiplication completed as float in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_print(m);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("sm"), "w");
            fprintf(f, "sm %.2f\n", a);
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%d\n", m->dimX);
            fprintf(f, "%d\n", m->dimY);
            for (int i = 0; i < m->dimY; i++) {
                for (int j = 0; j < m->dimX; j++) {
                    float v = matcoo_get(m, i, j);
                    fprintf(f, "%.2f ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_free(m);
    } else {
        matcoo_i *m = matcoo_i_fromfile(filepath);
        timer_start();
        m = matcoo_i_sm(m, a);
        double elapsed = timer_end();
        printf("Scalar multiplication completed as int in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_i_print(m);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("sm"), "w");
            fprintf(f, "sm %.2f\n", a);
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "int\n");
            fprintf(f, "%d\n", m->dimX);
            fprintf(f, "%d\n", m->dimY);
            for (int i = 0; i < m->dimY; i++) {
                for (int j = 0; j < m->dimX; j++) {
                    int v = matcoo_i_get(m, i, j);
                    fprintf(f, "%d ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_i_free(m);
    }
}

// runs the trace operation
void run_tr(char *filepath, bool should_print, bool should_log, int thread_count) {
    if (is_float(filepath)) {
        matcoo *m = matcoo_fromfile(filepath);
        timer_start();
        float t = matcoo_trace(m);
        double elapsed = timer_end();
        printf("Trace completed as float in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            printf("trace = %.2f\n", t);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("tr"), "w");
            fprintf(f, "tr\n");
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%.4f\n", t);
            fprintf(f, "%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_free(m);
    } else {
        matcoo_i *m = matcoo_i_fromfile(filepath);
        timer_start();
        int t = matcoo_i_trace(m);
        double elapsed = timer_end();
        printf("Trace completed as int in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            printf("trace = %d\n", t);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("tr"), "w");
            fprintf(f, "tr\n");
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "int\n");
            fprintf(f, "%d\n", t);
            fprintf(f, "%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_i_free(m);
    }
}

// runs the add operation
void run_ad(char *filepath1, char *filepath2, bool should_print, bool should_log, int thread_count) {
    if (is_float(filepath1)) {
        matcoo *m1 = matcoo_fromfile(filepath1);
        matcoo *m2 = matcoo_fromfile(filepath2);
        if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {
            printf("Error: Matrices must be of the same rank for operation addition.\n");
            exit(EXIT_FAILURE);
        }
        timer_start();
        matcoo *mres = matcoo_add(m1, m2);
        double elapsed = timer_end();
        printf("Add completed as float in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("ad"), "w");
            fprintf(f, "add\n");
            fprintf(f, "%s\n", filepath1);
            fprintf(f, "%s\n", filepath2);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    float v = matcoo_get(mres, i, j);
                    fprintf(f, "%.2f ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_free(m1);
        matcoo_free(m2);
    } else {
        matcoo_i *m1 = matcoo_i_fromfile(filepath1);
        matcoo_i *m2 = matcoo_i_fromfile(filepath2);
        if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {
            printf("Error: Matrices must be of the same rank for operation addition.\n");
            exit(EXIT_FAILURE);
        }
        timer_start();
        matcoo_i *mres = matcoo_i_add(m1, m2);
        double elapsed = timer_end();
        printf("Add completed as int in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_i_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("ad"), "w");
            fprintf(f, "add\n");
            fprintf(f, "%s\n", filepath1);
            fprintf(f, "%s\n", filepath2);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "int\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    int v = matcoo_i_get(mres, i, j);
                    fprintf(f, "%d ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_i_free(m1);
        matcoo_i_free(m2);
    }
}

// runs the transpose operation
void run_ts(char *filepath, bool should_print, bool should_log, int thread_count) {
    if (is_float(filepath)) {
        matcoo *m = matcoo_fromfile(filepath);
        timer_start();
        matcoo *mres = matcoo_transpose(m);
        double elapsed = timer_end();
        printf("Transpose completed as float in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("ts"), "w");
            fprintf(f, "ts\n");
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    float v = matcoo_get(mres, i, j);
                    fprintf(f, "%.2f ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_free(m);
    } else {
        matcoo_i *m = matcoo_i_fromfile(filepath);
        timer_start();
        matcoo_i *mres = matcoo_i_transpose(m);
        double elapsed = timer_end();
        printf("Transpose completed as int in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_i_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("ts"), "w");
            fprintf(f, "ts\n");
            fprintf(f, "%s\n", filepath);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    float v = matcoo_i_get(mres, i, j);
                    fprintf(f, "%.2f ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcoo_i_free(m);
    }
}

// runs the multiplication operation
void run_mm(char *filepath1, char *filepath2, bool should_print, bool should_log, int thread_count) {
    if (is_float(filepath1)) {
        matcsr *m1 = matcsr_fromfile(filepath1);
        matcsc *m2 = matcsc_fromfile(filepath2);

        if (m1->dimX != m2->dimY) {
            printf("Multiplication is only defined for two matrices where A.dimX == B.dimY\n");
            exit(EXIT_FAILURE);
        }

        timer_start();
        matcoo *mres = mat_multiply(m1, m2, thread_count);
        double elapsed = timer_end();
        printf("Multiplication completed as float in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("mm"), "w");
            fprintf(f, "mm\n");
            fprintf(f, "%s\n", filepath1);
            fprintf(f, "%s\n", filepath2);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "float\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    float v = matcoo_get(mres, i, j);
                    fprintf(f, "%.2f ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }

        matcsr_free(m1);
        matcsc_free(m2);
    } else {
        matcsr_i *m1 = matcsr_i_fromfile(filepath1);
        matcsc_i *m2 = matcsc_i_fromfile(filepath2);

        if (m1->dimX != m2->dimY) {
            printf("Multiplication is only defined for two matrices where A.dimX == B.dimY\n");
            exit(EXIT_FAILURE);
        }

        timer_start();
        matcoo_i *mres = mat_multiply_i(m1, m2, thread_count);
        double elapsed = timer_end();
        printf("Multiplication completed as int in %.4fs using %d threads.\n", elapsed, thread_count);
        if (should_print) {
            matcoo_i_print(mres);
        }

        if (should_log) {
            FILE *f = fopen(get_log_name("mm"), "w");
            fprintf(f, "mm\n");
            fprintf(f, "%s\n", filepath1);
            fprintf(f, "%s\n", filepath2);
            fprintf(f, "%d\n", thread_count);
            fprintf(f, "int\n");
            fprintf(f, "%d\n", mres->dimX);
            fprintf(f, "%d\n", mres->dimY);
            for (int i = 0; i < mres->dimY; i++) {
                for (int j = 0; j < mres->dimX; j++) {
                    int v = matcoo_i_get(mres, i, j);
                    fprintf(f, "%d ", v);
                }
            }
            fprintf(f, "\n%.4f\n", elapsed);
            fclose(f);
        }
    }
}

void print_usage() {
    printf("Usage:\n"
           "\tOperations:\n"
           "\t--sm \tScalar multiplication\n"
           "\t--tr \tTrace\n"
           "\t--ad \tAddition\n"
           "\t--ts \tTranspose\n"
           "\t--mm \tMatrix multiplication\n"
           "\n"
           "\tOptionals:\n"
           "\t-p  \tPrints result to stdout");
}

int main(int argc, char *argv[]) {

    //test_correctness("data/float123.in", "data/float123.in");

    // gather info about the arguments
    char operation = '0';
    int file_count = 1;
    char *file1 = NULL;
    char *file2 = NULL;
    float scalar;
    bool should_print = false;
    bool should_log = false;
    int thread_count = 8;
    // scalar multiplication
    if (argc < 2) {
        print_usage();
    } else {
        // go over every arg
        for (int i = 0; i < argc; i++) {
            char *v = argv[i];
            if (strcmp(v, "--sm") == 0) {
                operation = 's';
                scalar = strtod(argv[++i], NULL); // if the next arg is a file, then its not an option
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
                file1 = argv[++i]; // if the next arg is a file, then its not an option
                if (file_count > 1) {
                    file2 = argv[++i]; // handle a second file
                }
            } else if (strcmp(v, "-p") == 0) {
                should_print = true;
            } else if (strcmp(v, "-l") == 0) {
                should_log = true;
            } else if (strcmp(v, "-t") == 0) {
                thread_count = (int) strtol(argv[++i], NULL, 10); // if the next arg is a file, then its not an option
            }
        }
    }

    if (operation == '0') {
        printf("Error: Unknown operation\n");
        print_usage();
    } else if (!file1) {
        printf("Error: No files supplied, supply files with -f a.in [ b.in ]\n");
        print_usage();
    }

    // run an operation based on the arguments
    switch (operation) {
        case 's':
            run_sm(file1, scalar, should_print, should_log, thread_count);
            break;
        case 'r':
            run_tr(file1, should_print, should_log, thread_count);
            break;
        case 'a':
            run_ad(file1, file2, should_print, should_log, thread_count);
            break;
        case 't':
            run_ts(file1, should_print, should_log, thread_count);
            break;
        case 'm':
            run_mm(file1, file2, should_print, should_log, thread_count);
            break;
        default:
            print_usage();
            return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}