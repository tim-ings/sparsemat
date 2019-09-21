#ifndef SPARSEMAT_LIST_H
#define SPARSEMAT_LIST_H
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>


typedef struct list_ list;
struct list_ {
    float** parts;
    int parts_index; // the end of the parts array, used to append new parts
    int index; // the end of the array, used to append values
    int capacity; // the current max length of all the parts (i.e. the list)
    int increment; // the length of each part
};

list* list_new(int size, int partsize);
float list_get(list* l, int i);
void list_set(list* l, int i, float val);
void list_append(list** l, float val);

#endif
