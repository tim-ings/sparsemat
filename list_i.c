#include "list_i.h"


list_i *list_i_new(int capacity, int increment) {
    list_i *l = malloc(sizeof(list_i) * capacity);
    l->parts = malloc(sizeof(int*) * LIST_PART_COUNT);
    l->parts_length = 0;
    // init some initial parts so we meet requested capacity
    for (int i = 0; i < capacity / increment + 1; i++) {
        l->parts[i] = malloc(sizeof(int) * increment);
        l->parts_length++;
    }
    l->length = 0;
    l->capacity = capacity;
    l->increment = increment;
    return l;
}

void list_i_free(list_i *l) {
    for (int i = 0; i < l->increment; i++) { // each part is inc long
        if (l->parts[i])
            free(l->parts[i]);
    }
    free(l->parts);
    free(l);
}

int list_i_get(list_i *l, int i) {
    int pi = i / l->increment; // part index, int division floors result
    int si = i - (l->increment * pi); // sub index, the index within the part
    return l->parts[pi][si];
}

void list_i_set(list_i *lst, int i, int val) {
//    list_i* lst = *l;
    int pi = i / lst->increment; // part index, int division floors result
    int si = i - (lst->increment * pi); // sub index, the index within the part
    lst->parts[pi][si] = val;
}

void list_i_append(list_i **l, int val) {
    list_i *lst = *l;
    if (lst->length >= lst->increment * lst->parts_length) { // expand the array if we need to
        list_i_expand(lst);
    }
    list_i_set(lst, lst->length++, val);
}

void list_i_expand(list_i *lst) {
    lst->parts[lst->parts_length++] = (int *) malloc(sizeof(int) * lst->increment);
    lst->capacity += lst->increment;
}

void list_i_print(list_i *l) {
    for (int i = 0; i < l->length; i++) {
        printf("%d, ", list_i_get(l, i));
    }
    printf("\n");
}
