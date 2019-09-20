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
    matcoo* m = matcoo_blank();
    m->dimX = dimX;
    m->dimY = dimY;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float d = data[i * dimY + j];
            if (d != 0) {
                ll_float_push(m->coords_val, d * a);
                ll_float_push(m->coords_i, (float)i);
                ll_float_push(m->coords_j, (float)j);
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
    matcoo* m = matcoo_blank();
    m->dimX = dimX1;
    m->dimY = dimY1;
    for (int i = 0; i < dimY1; i++) {
        for (int j = 0; j < dimY1; j++) {
            int di = i * dimY1 + j;
            float d1 = data1[di];
            float d2 = data2[di];
            if (d1 != 0 || d2 != 0) {
                ll_float_push(m->coords_val, d1 + d2);
                ll_float_push(m->coords_i, (float)i);
                ll_float_push(m->coords_j, (float)j);
            }
        }
    }
    return m;
}

matcoo* transpose(const char* filepath) {
    int dimX;
    int dimY;
    float* data = readfile(filepath, &dimX, &dimY);
    matcoo* m = matcoo_blank();
    m->dimX = dimY; // swap x and y dims on the output mat
    m->dimY = dimX;
    for (int i = 0; i < dimY; i++) {
        for (int j = 0; j < dimX; j++) {
            float d = data[i * dimY + j];
            if (d != 0) {
                ll_float_push(m->coords_val, d);
                ll_float_push(m->coords_i, (float)j); // swap the i and j coords to transpose the value
                ll_float_push(m->coords_j, (float)i);
            }
        }
    }
    return m;
}

matcoo*  multiply(const char* filepath1, const char* filepath2) {
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

int main(int argc, char *argv[]) {

//    int len = 10000 * 1000;
//    timer_start();
//    for (int i = 0; i < len; i++) {
//        free(malloc(1));
//        free(malloc(1));
//    }
//    float t1 = timer_end();
//    printf("Combined loop time: %.4f\n", t1);
//    timer_start();
//    for (int i = 0; i < len; i++) {
//        free(malloc(1));
//    }
//    for (int i = 0; i < len; i++) {
//        free(malloc(1));
//    }
//    float t2 = timer_end();
//    printf("Split loop time: %.4f\n", t2);




    const char* f1 = "data/float123.in";
    const char* f2 = "data/float123.in";

    int dx, dy;
    float* d = readfile(f1, &dx, &dy);
    matcoo* m = matcoo_new(d, dx, dy);
    printf("Initial Matrix:\n");
    matcoo_print(m);

    matcoo* m_sm = sm(f1, 2);
    printf("Scalar Multiplication:\n");
    matcoo_print(m_sm);

    float t = trace(f1);
    printf("Trace: %.2f\n", t);

    matcoo* m_add = add(f1, f2);
    printf("Matrix Addition:\n");
    matcoo_print(m_add);

    matcoo* m_transpose = transpose(f1);
    printf("Matrix Transposition:\n");
    matcoo_print(m_transpose);

    matcoo* m_multiply = multiply(f1, f2);
//    printf("Matrix Multiplication:\n");
//    matcoo_print(m_multiply);

    return 0;
}