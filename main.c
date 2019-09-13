#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "matcoo.h"
#include "matcsr.h"
#include "matcsc.h"


#define LOAD_COUNT 1
#define ITT_COUNT 1500000
#define ITT_COUNT_BIG ITT_COUNT * 100

#define TEST_SM 0
#define TEST_TR 0
#define TEST_AD 1
#define TEST_TS 0
#define TEST_MM 0



struct timespec tstart={0,0}, tend={0,0};

void timer_start() {
    clock_gettime(CLOCK_MONOTONIC, &tstart);
}

float timer_end() {
    clock_gettime(CLOCK_MONOTONIC, &tend);
    return ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
}

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
            default: {
                printf("Unexpected line in %s, ignoring it\n", fileName);
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

    /*
     * LOAD MATRICES
     */
    matcoo* mcoo;
    timer_start();
    for (int i = 0; i < LOAD_COUNT; i++) {
        mcoo = matcoo_new(data, dimX, dimY);
    }
    float elapsed = timer_end();
    printf("Loaded %d COO matrices in %.5fs\n", LOAD_COUNT, elapsed);
    matcoo_print(mcoo);

    matcsr* mcsr;
    timer_start();
    for (int i = 0; i < LOAD_COUNT; i++) {
        mcsr = matcsr_new(data, dimX, dimY);
    }
    elapsed = timer_end();
    printf("Loaded %d CSR matrices in %.5fs\n", LOAD_COUNT, elapsed);
    matcsr_print(mcsr);

    matcsc* mcsc;
    timer_start();
    for (int i = 0; i < LOAD_COUNT; i++) {
        mcsc = matcsc_new(data, dimX, dimY);
    }
    elapsed = timer_end();
    printf("Loaded %d CSC matrices in %.5fs\n", LOAD_COUNT, elapsed);
    matcsc_print(mcsc);

    /*
     * CALCULATE TRACE
     */
    if (TEST_TR) {
        float sum;
        timer_start();
        for (int i = 0; i < ITT_COUNT_BIG; i++) {
            sum = matcoo_trace(mcoo);
        }
        elapsed = timer_end();
        printf("Calculated trace of %d COO matrices in %.5fs\n", ITT_COUNT_BIG, elapsed);
        printf("Trace COO: %.2f\n", sum);

        timer_start();
        for (int i = 0; i < ITT_COUNT_BIG; i++) {
            sum = matcsr_trace(mcsr);
        }
        elapsed = timer_end();
        printf("Calculated trace of %d CSR matrices in %.5fs\n", ITT_COUNT_BIG, elapsed);
        printf("Trace CSR: %.2f\n", sum);

        timer_start();
        for (int i = 0; i < ITT_COUNT_BIG; i++) {
            sum = matcsc_trace(mcsc);
        }
        elapsed = timer_end();
        printf("Calculated trace of %d CSC matrices in %.5fs\n", ITT_COUNT_BIG, elapsed);
        printf("Trace CSC: %.2f\n", sum);
    }

    /*
     * CALCULATE SCALAR MULTIPLE
     */
    if (TEST_SM) {
        int sm_count = (ITT_COUNT % 2 == 0) ? ((ITT_COUNT_BIG) + 1) : (ITT_COUNT_BIG);
        timer_start();
        for (int i = 0; i < sm_count; i++) {
            mcoo = matcoo_sm(mcoo, (i % 2 == 0) ? (10.0f) : (1.0f / 10.0f));
        }
        elapsed = timer_end();
        printf("Calculated scalar multiple of %d COO matrices in %.5fs\n", sm_count, elapsed);
        matcoo_print(mcoo);

        timer_start();
        for (int i = 0; i < sm_count; i++) {
            mcsr = matcsr_sm(mcsr, i % 2 == 0 ? 10.0f : 1.0f / 10.0f);
        }
        elapsed = timer_end();
        printf("Calculated scalar multiple of %d CSR matrices in %.5fs\n", sm_count, elapsed);
        matcsr_print(mcsr);

        timer_start();
        for (int i = 0; i < sm_count; i++) {
            mcsc = matcsc_sm(mcsc, i % 2 == 0 ? 10.0f : 1.0f / 10.0f);
        }
        elapsed = timer_end();
        printf("Calculated scalar multiple of %d CSC matrices in %.5fs\n", sm_count, elapsed);
        matcsr_print(mcsr);
    }

    if (TEST_AD) {
        printf("Adding m1 and m2:\n");
        matcoo_print(mcoo);
        printf("+\n");
        matcoo_print(mcoo);
        printf("=\n");
        matcoo* mad = matcoo_add(mcoo, mcoo);
        matcoo_print(mad);
    }

//    fclose(fp);
//    if (line)
//        free(line);

    return 0;
}

int main(int argc, char *argv[]) {
    readfile("data/float1.in");
    return 0;
}