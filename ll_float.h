#ifndef SPARSEMAT_LL_FLOAT_H
#define SPARSEMAT_LL_FLOAT_H
#include <stdlib.h>

#define PRINT_COUNT_MAX 15

typedef struct ll_float_ ll_float;
typedef struct ll_float_node_ ll_float_node;

struct ll_float_ {
    ll_float_node* first;
    ll_float_node* last;
    int length;
};

struct ll_float_node_ {
    float value;
    ll_float* list;
    ll_float_node* next;
    ll_float_node* prev;
};

void ll_float_push(ll_float* list, float value);
ll_float* ll_float_new();
ll_float_node* ll_float_get(ll_float* list, int index);
float ll_float_next(ll_float_node** cur);

#endif //SPARSEMAT_LL_FLOAT_H
