#include "list.h"


list* list_new(int capacity, int increment) {
    list* l = (list*)malloc(sizeof(list) * capacity);
    l->parts = (float**)malloc(sizeof(float*) * 100);
    // init some initial parts so we meet requested capacity
    for (int i = 0; i < capacity / increment + 1; i++) {
        l->parts[i] = (float*)malloc(sizeof(float) * increment);
    }
    l->length = 0;
    l->capacity = capacity;
    l->increment = increment;
    return l;
}

float list_get(list* l, int i) {
    int pi = i / l->increment; // part index, int division floors result
    int si = i - (l->increment * pi); // sub index, the index within the part
    return l->parts[pi][si];
}

void list_set(list* lst, int i, float val) {
//    list* lst = *l;
    int pi = i / lst->increment; // part index, int division floors result
    int si = i - (lst->increment * pi); // sub index, the index within the part
    lst->parts[pi][si] = val;
}

void list_append(list** l, float val) {
    list* lst = *l;
    if (lst->length >= lst->increment * lst->parts_index) { // expand the array if we need to
        list_expand(lst);
    }
    list_set(lst, lst->length++, val);
}

void list_expand(list* lst) {
    lst->parts[lst->parts_index++] = (float*)malloc(sizeof(float) * lst->increment);
    lst->capacity += lst->increment;
}

void list_print(list* l) {
    for (int i = 0; i < l->length; i++) {
        printf("%.2f, ", list_get(l, i));
    }
    printf("\n");
}

//void list_append(list* lst, float val) {
//    if (lst->length >= lst->increment * lst->parts_index) { // expand the array if we need to
//        list_expand(lst);
//    }
//    list_set(lst, lst->length++, val);
//}
