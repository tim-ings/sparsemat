#include "matcoo.h"


matcoo* matcoo_blank() {
    matcoo* m = malloc(sizeof(matcoo));
    m->coords_val = ll_float_new();
    m->coords_i = ll_float_new();
    m->coords_j = ll_float_new();
    m->dimX = 0;
    m->dimY = 0;
    return m;
}

void matcoo_free(matcoo* m) {
    ll_float_free(m->coords_val);
    ll_float_free(m->coords_i);
    ll_float_free(m->coords_j);
}

matcoo* matcoo_new(const float* data, int dimX, int dimY) {
    ll_float* coords_val = ll_float_new();
    ll_float* coords_i = ll_float_new();
    ll_float* coords_j = ll_float_new();
    if (data) {
        int stride = dimY;
        for (int i = 0; i < dimY; i++) {
            for (int j = 0; j < dimX; j++) {
                float val = data[i * stride + j];
                if (val != 0) {
                    ll_float_push(coords_val, val);
                    ll_float_push(coords_i, (float)i);
                    ll_float_push(coords_j, (float)j);
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

    return m;
}

void matcoo_print(matcoo* m) {
    for (int i = 0; i < m->dimY; i++) {
        for (int j = 0; j < m->dimX; j++) {
            float v = matcoo_get(m, i, j);
            printf("%.2f\t\t", v);
        }
        printf("\n");
    }
}

float matcoo_get(matcoo* m, int mi, int mj) {
    ll_float_node* cv = m->coords_val->first;
    ll_float_node* ci = m->coords_i->first;
    ll_float_node* cj = m->coords_j->first;
    while (cv) {
        if ((int)ci->value == mi && (int)cj->value == mj) {
            return cv->value;
        }
        cv = cv->next;
        ci = ci->next;
        cj = cj->next;
    }
    return 0.0f;
}
