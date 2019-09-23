#include "matcoo_i.h"

matcoo_i *matcoo_i_zeroes(int dx, int dy) {
    matcoo_i *m = malloc(sizeof(matcoo_i));
    m->vals = ll_int_new();
    m->is = ll_int_new();
    m->js = ll_int_new();
    m->dimX = dx;
    m->dimY = dy;
    return m;
}

void matcoo_i_free(matcoo_i *m) {
    ll_int_free(m->vals);
    ll_int_free(m->is);
    ll_int_free(m->js);
    free(m);
}

matcoo_i *matcoo_i_new(const int *data, int dimX, int dimY) {
    // init a new matcoo_i struct
    matcoo_i *m = (matcoo_i *) malloc(sizeof(matcoo_i));
    m->dimX = dimX;
    m->dimY = dimY;
    // each non zro value is stored as a (val, i, j)
    m->vals = ll_int_new();
    m->is = ll_int_new();
    m->js = ll_int_new();
    if (data) {
        int stride = dimY;
        for (int i = 0; i < dimY; i++) {
            for (int j = 0; j < dimX; j++) {
                int val = data[i * stride + j];
                if (val != 0) {
                    ll_int_push(m->vals, val);
                    ll_int_push(m->is, i);
                    ll_int_push(m->js, j);
                }
            }
        }
        ll_int_push(m->vals, data[dimY * stride + dimX]);
        ll_int_push(m->is, dimY);
        ll_int_push(m->js, dimX);
    }
    return m;
}

matcoo_i *matcoo_i_build(matcoo_i *m, int val, int i, int j) {
    if (val != 0) {
        ll_int_push(m->vals, val);
        ll_int_push(m->is, i);
        ll_int_push(m->js, j);
    }
    return m;
}

void matcoo_i_print(matcoo_i *m) {
    for (int i = 0; i < m->dimY; i++) {
        for (int j = 0; j < m->dimX; j++) {
            int v = matcoo_i_get(m, i, j);
            printf("%d\t\t", v);
        }
        printf("\n");
    }
}

int matcoo_i_get(matcoo_i *m, int i, int j) {
    ll_int_node *cv = m->vals->first;
    ll_int_node *ci = m->is->first;
    ll_int_node *cj = m->js->first;
    while (cv) {
        if ((int) ci->value == i && (int) cj->value == j) {
            return cv->value;
        }
        cv = cv->next;
        ci = ci->next;
        cj = cj->next;
    }
    return 0.0f;
}

matcoo_i *matcoo_i_sm(matcoo_i *m, float a) {
    ll_int_node *c = m->vals->first;
    while (c) {
        c->value *= a;
        c = c->next;
    }
    return m;
}

int matcoo_i_trace(matcoo_i *m) {
    if (m->dimY != m->dimX) {
        assert("Trace is only well defined for square matrices." && 0);
    }
    ll_int_node *cv = m->vals->first;
    ll_int_node *ci = m->is->first;
    ll_int_node *cj = m->js->first;
    int sum = 0;
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

matcoo_i *matcoo_i_add(matcoo_i *m1, matcoo_i *m2) {
    ll_int_node *c1 = m1->vals->first;
    ll_int_node *i1 = m1->is->first;
    ll_int_node *j1 = m1->js->first;
    ll_int_node *c2 = m2->vals->first;
    ll_int_node *i2 = m2->is->first;
    ll_int_node *j2 = m2->js->first;
    // go over every non-0 value in m1
    while (c1) {
        // go over every non-0 value in m2
        while (c2) {
            // if we have a match, add them and store in m1
            if (i1->value == i2->value && j1->value == j2->value) {
                c1->value += c2->value;
                // remove the used value from m2 so we dont check them again
                ll_int_remove(c2);
                ll_int_remove(i2);
                ll_int_remove(j2);
                // reset to start of m2 and break
                c2 = m2->vals->first;
                i2 = m2->is->first;
                j2 = m2->js->first;
                break;
            }
            c2 = c2->next;
            i2 = i2->next;
            j2 = j2->next;
        }
        c1 = c1->next;
        i1 = i1->next;
        j1 = j1->next;
    }
    return m1;
}

matcoo_i *matcoo_i_transpose(matcoo_i *m) {
    ll_int_node *i = m->is->first;
    ll_int_node *j = m->js->first;
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

matcoo_i *matcoo_i_multiply(matcoo_i *m1, matcoo_i *m2) {
    if (m1->dimY != m2->dimX) {
        assert("Multiplication is only defined for matrices where first.dimY == second.dimX" && 0);
    }
    matcoo_i *res = matcoo_i_zeroes(m2->dimX, m1->dimY);
    for (int i = 0; i < m1->dimY; i++) {
        for (int j = 0; j < m2->dimX; j++) {
            float sum = 0;
            for (int k = 0; k < m1->dimX; k++) {
                int v1 = matcoo_i_get(m1, i, k);
                int v2 = matcoo_i_get(m2, k, j);
                sum += v1 * v2;
            }
            matcoo_i_build(res, sum, i, j);
        }
    }
    return res;
}

bool matcoo_i_equals(matcoo_i *m1, matcoo_i *m2) {
    if (m1->dimX != m2->dimX || m1->dimY != m2->dimY) {
        return false;
    }
    ll_int_node *v1 = m1->vals->first;
    ll_int_node *i1 = m1->is->first;
    ll_int_node *j1 = m1->js->first;
    ll_int_node *v2 = m2->vals->first;
    ll_int_node *i2 = m2->is->first;
    ll_int_node *j2 = m2->js->first;
    while (v1 && v2) {
        if (v1->value != v2->value) {
            return false;
        }
        v1 = v1->next;
        i1 = i1->next;
        j1 = j1->next;
        v2 = v2->next;
        i2 = i2->next;
        j2 = j2->next;
    }
    return true;
}
