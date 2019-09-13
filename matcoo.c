#include "matcoo.h"


matcoo* matcoo_new(const float* data, int dimX, int dimY) {
#ifdef DEBUG_VERBOSE
    printf("\nBuilding new %dx%d COO matrix\n[ ", dimX, dimY);
#endif
    ll_float* coords_val = ll_float_new();
    ll_float* coords_i = ll_float_new();
    ll_float* coords_j = ll_float_new();
    if (data) {
        int stride = dimY;
#ifdef DEBUG_VERBOSE
        int prints_so_far = 0;
#endif
        for (int i = 0; i < dimY; i++) {
            for (int j = 0; j < dimX; j++) {
                float val = data[i * stride + j];
                if (val != 0) {
                    ll_float_push(coords_val, val);
                    ll_float_push(coords_i, (float)i);
                    ll_float_push(coords_j, (float)j);
#ifdef DEBUG_VERBOSE
                    if (prints_so_far < PRINT_COUNT_MAX / 1.5) {
                    printf("(%d, %d, %.2f), ", j, i, val);
                    prints_so_far++;
                }
#endif
                }
            }
        }
        ll_float_push(coords_val, data[dimY * stride + dimX]);
        ll_float_push(coords_i, (float)dimY);
        ll_float_push(coords_j, (float)dimX);
    }

    matcoo* m = (matcoo*)malloc(sizeof(matcoo));
    m->dimX = dimX;
    m->dimY = dimY;
    m->coords_val = coords_val;
    m->coords_i = coords_i;
    m->coords_j = coords_j;

#ifdef DEBUG_VERBOSE
    printf("]...\n%dx%d COO matrix built\n", dimX, dimY);
#endif
    return m;
}

matcoo* matcoo_sm(matcoo* m, float s) {
    ll_float_node* cur_val = m->coords_val->first;
    while (cur_val->next != NULL) {
        cur_val->value *= s;
        cur_val = cur_val->next;
    }
    return m;
}

void matcoo_print(matcoo* m) {
    ll_float_node* cur_val = m->coords_val->first;
    ll_float_node* cur_i = m->coords_i->first;
    ll_float_node* cur_j = m->coords_j->first;
    for (int i = 0; i < m->dimY; i++) {
        for (int j = 0; j < m->dimY; j++) {
            if (i == (int)cur_i->value && j == (int)cur_j->value) {
                printf("%.2f\t\t", cur_val->value);
                cur_val = cur_val->next;
                cur_i = cur_i->next;
                cur_j = cur_j->next;
            } else {
                printf("0.00\t\t");
            }
        }
        printf("\n");
    }
}

float matcoo_trace(matcoo* m) {
    if (m->dimY == m->dimX) {
        float sum = 0;
        ll_float_node* cur_val = m->coords_val->first;
        ll_float_node* cur_i = m->coords_i->first;
        ll_float_node* cur_j = m->coords_j->first;
        for (int i = 0; i < m->coords_val->length; i++) {
            if ((int)cur_i->value == (int)cur_j->value) {
                sum += cur_val->value;
            }
            cur_val = cur_val->next;
            cur_i = cur_i->next;
            cur_j = cur_j->next;
        }
        return sum;
    } else {
        printf("Trace is only well defined for square matrices.\n");
        return -1.0f;
    }
}

matcoo* matcoo_add(matcoo* m1, matcoo* m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {
        printf("Matrices must have the same rank (%d->%d, %d->%d)\n", m1->dimX, m2->dimX, m1->dimY, m2->dimY);
        return NULL;
    }
    matcoo* nm = matcoo_new(NULL, m1->dimX, m1->dimY);
    for (int i = 0; i < m1->dimX; i++) {
        for (int j = 0; j < m1->dimY; j++) {
            ll_float_node* n1 = matcoo_get(m1, i, j);
            ll_float_node* n2 = matcoo_get(m2, i, j);
            if (n1 && n2) {
                matcoo_set(nm, i, j, n1->value + n2->value);
            }
        }
    }
    return nm;
}

ll_float_node* matcoo_get(matcoo* m, int i, int j) {
    ll_float_node* cur_val = m->coords_val->first;
    ll_float_node* cur_i = m->coords_i->first;
    ll_float_node* cur_j = m->coords_j->first;
    for (int k = 0; k < m->coords_val->length; k++) {
        if (cur_i && (int)cur_i->value == i && cur_j && (int)cur_j->value == j) {
            return cur_val;
        }
        cur_val = cur_val->next;
        cur_i = cur_i->next;
        cur_j = cur_j->next;
    }
    return NULL;
}

void matcoo_set(matcoo* m, int i, int j, float val) {
    ll_float_node* cur_val = m->coords_val->first;
    ll_float_node* cur_i = m->coords_i->first;
    ll_float_node* cur_j = m->coords_j->first;
    if (m->coords_val->length > 0) {
        for (int k = 0; k < m->coords_val->length; k++) {
            if (cur_i && (int)cur_i->value == i && cur_j && (int)cur_j->value == j) {
                cur_val->value = val;
                return;
            }
            cur_val = cur_val->next;
            cur_i = cur_i->next;
            cur_j = cur_j->next;
        }
    }
    ll_float_push(m->coords_val, val);
    ll_float_push(m->coords_i, (float)i);
    ll_float_push(m->coords_j, (float)j);
}
