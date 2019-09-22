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
            printf("checking i,j pair (%d, %d) == (%d, %d)\n", (int)i1->value, (int)j1->value, (int)i2->value, (int)j2->value);
            // if we have a match, add them and store in m1
            if (i1->value == i2->value && j1->value == j2->value) {
                printf("\tfound a matching i,j pair (%d, %d)\n", (int)i1->value, (int)j1->value);
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
