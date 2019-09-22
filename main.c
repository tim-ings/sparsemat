#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "matcoo.h"
#include "matcsr.h"
#include "matcsc.h"


struct timespec tstart={0,0}, tend={0,0};

void timer_start() {
    clock_gettime(CLOCK_MONOTONIC, &tstart);
}

float timer_end() {
    clock_gettime(CLOCK_MONOTONIC, &tend);
    return ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
}

float* readfile(const char* fileName, int* out_dimX, int* out_dimY) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s\n", fileName);
        return NULL;
    }
    int counter = 0;
    float* data = NULL;
    while ((read = getline(&line, &len, fp)) != -1) {
        line[read - 1] = '\0';
        switch (counter) {
            case 0: {
                break;
            }
            case 1: {
                *out_dimX = (int)strtol(line, NULL, 10);
                break;
            }
            case 2: {
                *out_dimY = (int)strtol(line, NULL, 10);
                data = (float*)malloc(sizeof(float) * (*out_dimX) * (*out_dimY));
                break;
            }
            case 3: {
                const char delim[2] = " ";
                char* token = strtok(line, delim);
                char* err;
                data[0] = (float)strtod(token, &err);
                for (int i = 1; i < ((*out_dimX) * (*out_dimY)); i++) {
                    token = strtok(NULL, delim);
                    data[i] = (float)strtod(token, &err);
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

matcoo* matcoo_fromfile(const char* fp) {
    int dx, dy;
    float* d = readfile(fp, &dx, &dy);
    matcoo* m = matcoo_zeroes(dx, dy);

    for (int i = 0; i < dy; i++) {
        for (int j = 0; j < dx; j++) {
            float v = d[i * dx + j];
            m = matcoo_build(m, v, i, j);
        }
    }
    return m;
}

int main(int argc, char *argv[]) {
    const char* f1 = "data/float1.in";
    const char* f2 = "data/float1.in";

    printf("\tinit:\n");
    matcoo* m_init = matcoo_fromfile(f1);
    matcoo_print(m_init);
    matcoo_free(m_init);

    printf("\tsm:\n");
    matcoo* m_sm = matcoo_fromfile(f1);
    matcoo* m_sm_res = matcoo_sm(m_sm, 2);
    matcoo_print(m_sm_res);
    matcoo_free(m_sm_res);

    printf("\ttrace:\n");
    matcoo* m_trace = matcoo_fromfile(f1);
    float t = matcoo_trace(m_trace);
    printf("t=%.2f\n", t);
    matcoo_free(m_trace);

    printf("\tadd:\n");
    matcoo* m_add1 = matcoo_fromfile(f1);
    matcoo* m_add2 = matcoo_fromfile(f2);
    matcoo* m_add_res = matcoo_add(m_add1, m_add2);
    matcoo_print(m_add_res);
    matcoo_free(m_add_res);
    matcoo_free(m_add2);

    printf("\ttranspose:\n");
    matcoo* m_transpose = matcoo_fromfile(f1);
    matcoo* m_transpose_res = matcoo_transpose(m_transpose);
    matcoo_print(m_transpose_res);
    matcoo_free(m_transpose_res);

    printf("\tmultiply:\n");
    matcoo* m_multiply1 = matcoo_fromfile("data/float_dp23.in");
    matcoo* m_multiply2 = matcoo_fromfile("data/float_dp34.in");
    printf("m_multiply1:\n");
    matcoo_print(m_multiply1);
    printf("m_multiply2:\n");
    matcoo_print(m_multiply2);
    matcoo* m_multiply_res = matcoo_multiply(m_multiply1, m_multiply2);
    printf("m_multiply_res:\n");
    matcoo_print(m_multiply_res);
    matcoo_free(m_multiply_res);

    return 0;
}