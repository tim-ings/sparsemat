#include "ll_float.h"


void ll_float_push(ll_float* list, float value) {
    ll_float_node* new_node = malloc(sizeof(ll_float_node));
    new_node->list = list;
    new_node->next = NULL;
    new_node->value = value;
    if (list->length > 0) {
        new_node->prev = list->last;
        list->last->next = new_node;
        list->last = new_node;
    } else {
        list->last = new_node;
        list->first = new_node;
        new_node->prev = NULL;
    }
    list->length++;
}

ll_float* ll_float_new() {
    ll_float* list = (ll_float*)malloc(sizeof(ll_float));
    list->length = 0;
    list->first = NULL;
    list->last = NULL;
    return list;
}

ll_float_node* ll_float_get(ll_float* list, int index) {
    ll_float_node* current = list->first;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    return current;
}

float ll_float_next(ll_float_node** cur) {
    float val = (*cur)->value;
    *cur = (*cur)->next;
    return val;
}
