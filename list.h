#ifndef SPARSEMAT_LIST_H
#define SPARSEMAT_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define LIST_PART_COUNT 100

typedef struct list_ list;
struct list_ {
    float **parts;
    int parts_length;
    int length; // the current length of the list (sum of all parts)
    int capacity; // the current max length of the list (sum of all parts)
    int increment; // the length of each part
};

list *list_new(int capacity, int increment);

void list_free(list* l);

float list_get(list *l, int i);

void list_set(list *l, int i, float val);

void list_append(list **l, float val);

void list_expand(list *l);

void list_print(list *l);

#endif
