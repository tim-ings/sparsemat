#include "matcoo.h"


matcoo* matcoo_zeroes(int dx, int dy) {
    matcoo* m = malloc(sizeof(matcoo));
    m->vals = ll_float_new();
    m->is = ll_float_new();
    m->js = ll_float_new();
    m->dimX = dx;
    m->dimY = dy;
    return m;
}

void matcoo_free(matcoo* m) {
    ll_float_free(m->vals);
    ll_float_free(m->is);
    ll_float_free(m->js);
    free(m);
}

matcoo* matcoo_new(const float* data, int dimX, int dimY) {
    // init a new matcoo struct
    matcoo* m = (matcoo*)malloc(sizeof(matcoo));
    m->dimX = dimX;
    m->dimY = dimY;
    // each non zro value is stored as a (val, i, j)
    m->vals = ll_float_new();
    m->is = ll_float_new();
    m->js = ll_float_new();
    if (data) {
        int stride = dimY;
        for (int i = 0; i < dimY; i++) {
            for (int j = 0; j < dimX; j++) {
                float val = data[i * stride + j];
                if (val != 0) {
                    ll_float_push(m->vals, val);
                    ll_float_push(m->is, (float)i);
                    ll_float_push(m->js, (float)j);
                }
            }
        }
        ll_float_push(m->vals, data[dimY * stride + dimX]);
        ll_float_push(m->is, (float)dimY);
        ll_float_push(m->js, (float)dimX);
    }
    return m;
}

matcoo* matcoo_build(matcoo* m, float val, int i, int j) {
    if (val != 0) {
        ll_float_push(m->vals, val);
        ll_float_push(m->is, (float)i);
        ll_float_push(m->js, (float)j);
    }
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
    ll_float_node* cv = m->vals->first;
    ll_float_node* ci = m->is->first;
    ll_float_node* cj = m->js->first;
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
