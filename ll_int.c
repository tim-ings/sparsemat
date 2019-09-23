#include "ll_int.h"


void ll_int_push(ll_int *list, int value) {
    ll_int_node *new_node = malloc(sizeof(ll_int_node));
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

ll_int *ll_int_new() {
    ll_int *list = (ll_int *) malloc(sizeof(ll_int));
    list->length = 0;
    list->first = NULL;
    list->last = NULL;
    return list;
}

ll_int_node *ll_int_get(ll_int *list, int index) {
    ll_int_node *current = list->first;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    return current;
}

int ll_int_next(ll_int_node **cur) {
    int val = (*cur)->value;
    *cur = (*cur)->next;
    return val;
}

void ll_int_free(ll_int *ll) {
    ll_int_node *c = ll->first;
    while (c) {
        ll_int_node *old = c;
        ll_int_next(&c);
        free(old);
    }
    free(ll);
}

void ll_int_remove(ll_int_node *n) {
    if (n->prev) {
        n->prev->next = n->next;
    } else {
        n->list->first = n->next;
    }
    if (n->next) {
        n->next->prev = n->prev;
    } else {
        n->list->last = n->prev;
    }
}
