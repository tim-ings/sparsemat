#include "list.h"


list* list_new(int size, int partsize) {
    list* l = (list*)malloc(sizeof(list) * size);
    l->parts = (float**)malloc(sizeof(float*) * 100);
    l->index = 0;
    l->capacity = size;
    l->increment = partsize;
    return l;
}

float list_get(list* l, int i) {
    int pi = i / l->increment; // part index, int division floors result
    int si = i - (l->increment * pi); // sub index, the index within the part
    return l->parts[pi][si];
}

void list_set(list* l, int i, float val) {
    int pi = i / l->increment; // part index, int division floors result
    int si = i - (l->increment * pi); // sub index, the index within the part
    l->parts[pi][si] = val;
}

void list_append(list** l, float val) {
    list* lst = *l;
    if (lst->index >= lst->increment * lst->parts_index) { // expand the array if we need to
        lst->parts[lst->parts_index++] = (float*)malloc(sizeof(float) * lst->increment);
        lst->capacity += lst->increment;
    }
    list_set(*l, lst->index++, val);

}
