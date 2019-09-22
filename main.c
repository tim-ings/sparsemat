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

float** readfile_arr(const char* fileName, int* out_dimX, int* out_dimY) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s", fileName);
        return NULL;
    }
    int counter = 0;
    float** data = NULL;
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
                data = (float**)malloc(sizeof(float*) * (*out_dimY));
                break;
            }
            case 3: {
                char* err;
                char* token;
                const char delim[2] = " ";
                for (int i = 0; i < *out_dimY; i++) {
                    data[i] = malloc(sizeof(float) * (*out_dimX));
                    for (int j = 0; j < *out_dimX; j++) {
                        token = strtok(i == 0 && j == 0 ? line : NULL, delim);
                        data[i][j] = (float)strtof(token, &err);
                    }
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

float* readfile(const char* fileName, int* out_dimX, int* out_dimY) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s", fileName);
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

matcoo* sm(const char* filepath, float a) {
    int dimX;
    int dimY;
    float* data = readfile(filepath, &dimX, &dimY);
    matcoo* m = matcoo_zeroes(dimX, dimY);
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float d = data[i * dimY + j];
            if (d != 0) {
                ll_float_push(m->vals, d * a);
                ll_float_push(m->is, (float)i);
                ll_float_push(m->js, (float)j);
            }
        }
    }
    return m;
}

float trace(const char* filepath) {
    int dimX;
    int dimY;
    float* data = readfile(filepath, &dimX, &dimY);
    float diagsum = 0;
    for (int i = 0; i < dimX; i++) {
        diagsum += data[i * (dimX + 1)];
    }
    return diagsum;
}

matcoo* add(const char* filepath1, const char* filepath2) {
    int dimX1;
    int dimY1;
    float* data1 = readfile(filepath1, &dimX1, &dimY1);
    int dimX2;
    int dimY2;
    float* data2 = readfile(filepath2, &dimX2, &dimY2);
    if (dimX1 != dimX2 || dimY1 != dimY2) {
        return NULL;
    }
    matcoo* m = matcoo_zeroes(dimX1, dimY1);
    for (int i = 0; i < dimY1; i++) {
        for (int j = 0; j < dimY1; j++) {
            int di = i * dimY1 + j;
            float d1 = data1[di];
            float d2 = data2[di];
            if (d1 != 0 || d2 != 0) {
                ll_float_push(m->vals, d1 + d2);
                ll_float_push(m->is, (float)i);
                ll_float_push(m->js, (float)j);
            }
        }
    }
    return m;
}

matcoo* transpose(const char* filepath) {
    int dimX;
    int dimY;
    float* data = readfile(filepath, &dimX, &dimY);
    matcoo* m = matcoo_zeroes(dimX, dimY);
    m->dimX = dimY; // swap x and y dims on the output mat
    m->dimY = dimX;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float d = data[i * dimY + j];
            if (d != 0) {
                ll_float_push(m->vals, d);
                ll_float_push(m->is, (float)j); // swap the i and j coords to transpose the value
                ll_float_push(m->js, (float)i);
            }
        }
    }
    return m;
}

float**  multiply(const char* filepath1, const char* filepath2) {
    int dimX1;
    int dimY1;
    int dimX2;
    int dimY2;
    float** data1 = readfile_arr(filepath1, &dimX1, &dimY1);
    float** data2 = readfile_arr(filepath2, &dimX2, &dimY2);
    float** mul = (float**)malloc(sizeof(float*) * dimY1);
    for (int i = 0; i < dimY1; i++)
    {
        mul[i] = (float*)malloc(sizeof(float) * dimX1);
        for (int j = 0; j < dimX1; j++)
        {
            for (int k = 0; k < dimY1; k++)
            {
                mul[i][j] += data1[i][k] * data2[k][j];
            }
        }
    }

    printf("Matrix Multiplication:\n");
    for (int i = 0; i < dimY1; i++){
        for (int j = 0; j < dimX1; j++) {
            printf("%.2f\t\t", mul[i][j]);
        }
        printf("\n");
    }

    return mul;
}

matcoo* get_matcoo(const char* fp) {
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
    matcoo* m_init = get_matcoo(f1);
    matcoo_print(m_init);

    printf("\tsm:\n");
    matcoo* m_sm = get_matcoo(f1);
    matcoo* m_sm_res = matcoo_sm(m_sm, 2);
    matcoo_print(m_sm_res);

    printf("\ttrace:\n");
    matcoo* m_trace = get_matcoo(f1);
    float t = matcoo_trace(m_trace);
    printf("t=%.2f\n", t);

    printf("\tadd:\n");
    matcoo* m_add1 = get_matcoo(f1);
    matcoo* m_add2 = get_matcoo(f2);
    matcoo* m_add_res = matcoo_add(m_add1, m_add2);
    matcoo_print(m_add_res);

    printf("\ttranspose:\n");
    matcoo* m_transpose = get_matcoo(f1);
    matcoo* m_transpose_res = matcoo_transpose(m_transpose);
    matcoo_print(m_transpose_res);

    return 0;
}