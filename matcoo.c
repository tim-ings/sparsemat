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
        _matcoo_itt_next(&cv, &ci, &cj);
    }
    return 0.0f;
}

matcoo* matcoo_sm(matcoo* m, float a) {
    ll_float_node *c = m->vals->first;
    while (c) {
        c->value *= a;
        c = c->next;
    }
    return m;
}

float matcoo_trace(matcoo* m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    ll_float_node *cv = m->vals->first;
    ll_float_node *ci = m->is->first;
    ll_float_node *cj = m->js->first;
    float sum = 0.0f;
    while (cv) {
        if (ci->value == cj->value) {
            sum += cv->value;
        }
        cv = cv->next;
        ci = ci->next;
        cj = cj->next;
    }
    return sum;
}

matcoo* matcoo_add(matcoo* m1, matcoo* m2) {
    ll_float_node *c1 = m1->vals->first;
    ll_float_node *i1 = m1->is->first;
    ll_float_node *j1 = m1->js->first;
    ll_float_node *c2 = m2->vals->first;
    ll_float_node *i2 = m2->is->first;
    ll_float_node *j2 = m2->js->first;
    // go over every non-0 value in m1
    while (c1) {
        // go over every non-0 value in m2
        while (c2) {
            // if we have a match, add them and store in m1
            if (i1->value == i2->value && j1->value == j2->value) {
                c1->value += c2->value;
                // remove the used value from m2 so we dont check them again
                ll_float_remove(c2);
                ll_float_remove(i2);
                ll_float_remove(j2);
                // reset to start of m2 and break
                c2 = m2->vals->first;
                i2 = m2->is->first;
                j2 = m2->js->first;
                break;
            }
            _matcoo_itt_next(&c2, &i2, &j2);
        }
        _matcoo_itt_next(&c1, &i1, &j1);
    }
    return m1;
}

matcoo* matcoo_transpose(matcoo* m) {
    ll_float_node *i = m->is->first;
    ll_float_node *j = m->js->first;
    while (i) {
        // swap i and j vals to transpose
        float oldi = i->value;
        i->value = j->value;
        j->value = oldi;
        i = i->next;
        j = j->next;
    }
    return m;
}

void _matcoo_itt_next(ll_float_node** c, ll_float_node** i, ll_float_node** j) {
    *c = (*c)->next;
    *i = (*i)->next;
    *j = (*j)->next;
}

matcoo* matcoo_multiply(matcoo* m1, matcoo* m2) {
    if (m1->dimY != m2->dimX) {
        assert("Multiplication is only defined for matrices where first.dimY == second.dimX" && 0);
    }
    matcoo* res = matcoo_zeroes(m2->dimX, m1->dimY);
    for (int i = 0; i < m1->dimY; i++) {
        for (int j = 0; j < m2->dimX; j++) {
            float sum = 0;
            for (int k = 0; k < m1->dimX; k++) {
                float v1 = matcoo_get(m1, i, k);
                float v2 = matcoo_get(m2, k, j);
                sum += v1 * v2;
            }
            matcoo_build(res, sum, i, j);
        }
    }
    return res;
}
