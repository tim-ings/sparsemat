#ifndef SPARSEMAT_LIST_I_H
#define SPARSEMAT_LIST_I_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define LIST_PART_COUNT 100

typedef struct list_i_ list_i;
struct list_i_ {
    int **parts;
    int parts_length;
    int length; // the current length of the list (sum of all parts)
    int capacity; // the current max length of the list (sum of all parts)
    int increment; // the length of each part
};

list_i *list_i_new(int capacity, int increment);

void list_i_free(list_i* l);

int list_i_get(list_i *l, int i);

void list_i_set(list_i *l, int i, int val);

void list_i_append(list_i **l, int val);

void list_i_expand(list_i *l);

void list_i_print(list_i *l);

#endif
