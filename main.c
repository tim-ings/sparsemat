#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ll_float.h"
#include "matcoo.h"
#include "matcsr.h"
#include "matcsc.h"


int readfile(const char* fileName) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s", fileName);
        return 1;
    }
    int dimX = -1;
    int dimY = -1;
    float* data = NULL;
    int counter = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
#ifdef DEBUG_VERBOSE
        printf("Retrieved line of length %zu\n", read);
#endif
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
                data = (float*)malloc(sizeof(float) * dimX * dimY);
                break;
            }
            case 3: {
                const char delim[2] = " ";
                char* token = strtok(line, delim);
                char* err;
                data[0] = (float)strtod(token, &err);
                for (int i = 1; i < dimX * dimY; i++) {
                    token = strtok(NULL, delim);
                    data[i] = (float)strtod(token, &err);
                }
                break;
            }
        }
        counter++;
    }

#ifdef DEBUG_VERBOSE
    printf("DIM X = '%d'\n", dimX);
    printf("DIM Y = '%d'\n", dimY);

    for (int i = 0; i < dimX * dimY; i++) {
        printf("%f, ", data[i]);
    }
    printf("\n");
#endif
    printf("%f\n", data[0]);

    fclose(fp);
    if (line)
        free(line);

    printf("\nNEW COO:\n");
    matcoo* matcoo = matcoo_new(data, dimX, dimY);
    matcoo_print(matcoo);
    printf("\nNEW CSR:\n");
    matcsr* matcsr = matcsr_new(data, dimX, dimY);
    matcsr_print(matcsr);
    printf("\nNEW CSC:\n");
    matcsc* matcsc = matcsc_new(data, dimX, dimY);
    matcsc_print(matcsc);

    float sum = matcoo_trace(matcoo);
    printf("\nTrace COO: %.2f\n", sum);

    sum = matcsr_trace(matcsr);
    printf("\nTrace CSR: %.2f\n", sum);

    sum = matcsc_trace(matcsc);
    printf("\nTrace CSC: %.2f\n", sum);

    printf("\nSM COO * 10.0f:\n");
    matcoo = matcoo_sm(matcoo, 10.0f);
    matcoo_print(matcoo);

    printf("\nSM CSR * 10.0f:\n");
    matcsr = matcsr_sm(matcsr, 10.0f);
    matcsr_print(matcsr);

    printf("\nSM CSC * 10.0f:\n");
    matcsc = matcsc_sm(matcsc, 10.0f);
    matcsc_print(matcsc);

    return 0;
}

int main(int argc, char *argv[]) {
    readfile("data/float1.in");
    return 0;
}