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
        printf("Retrieved line of length %zu\n", read);
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

    printf("DIM X = '%d'\n", dimX);
    printf("DIM Y = '%d'\n", dimY);

    for (int i = 0; i < dimX * dimY; i++) {
        printf("%f, ", data[i]);
    }
    printf("\n");
    printf("\n");

    fclose(fp);
    if (line)
        free(line);

    matcoo* matcoo = matcoo_new(data, dimX, dimY);
    matcsr* matcsr = matcsr_new(data, dimX, dimY);
    matcsc* matcsc = matcsc_new(data, dimX, dimY);
    return 0;
}

int main(int argc, char *argv[]) {
    readfile("data/float64.in");
}