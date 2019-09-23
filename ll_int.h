#ifndef SPARSEMAT_LL_INT_H
#define SPARSEMAT_LL_INT_H

#include <stdlib.h>


typedef struct ll_int_ ll_int;
typedef struct ll_int_node_ ll_int_node;

struct ll_int_ {
    ll_int_node *first;
    ll_int_node *last;
    int length;
};

struct ll_int_node_ {
    int value;
    ll_int *list;
    ll_int_node *next;
    ll_int_node *prev;
};

void ll_int_push(ll_int *list, int value);

ll_int *ll_int_new();

int ll_int_next(ll_int_node **cur);

void ll_int_free(ll_int *ll);

void ll_int_remove(ll_int_node *n);

#endif
