#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int readfile(const char* fileName) {
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) {
        printf("File does not exist: %s", fileName);
        return 1;
    }
    float dimX = -1;
    float dimY = -1;
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
                dimX = (float)strtod(line, NULL);
                break;
            }
            case 2: {
                dimY = (float)strtod(line, NULL);
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

    printf("DIM X = '%f'\n", dimX);
    printf("DIM Y = '%f'\n", dimY);

    for (int i = 0; i < dimX * dimY; i++) {
        printf("%f, ", data[i]);
    }

    fclose(fp);
    if (line)
        free(line);

    return 0;
}

int main(int argc, char *argv[])
{
//    bool isCaseInsensitive = false;
//    int opt;
//    enum { CHARACTER_MODE, WORD_MODE, LINE_MODE } mode = CHARACTER_MODE;
//
//    while ((opt = getopt(argc, argv, "ilw")) != -1) {
//        switch (opt) {
//            case 'i': isCaseInsensitive = true; break;
//            case 'l': mode = LINE_MODE; break;
//            case 'w': mode = WORD_MODE; break;
//            default:
//                fprintf(stderr, "Usage: %s [-ilw] [file...]\n", argv[0]);
//                exit(EXIT_FAILURE);
//        }
//    }

    readfile("../data/float1.in");
}