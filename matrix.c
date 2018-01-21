#include "matrix.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _mat_get(mat, i, j) ((mat)->v[(i) * ((mat)->max_m) + (j)])
#define _mat_set(mat, i, j, val) ((mat)->v[(i) * ((mat)->max_m) + (j)] = (val))
// #define _mat_v_add(mat,i,j,a) ((mat)->v[(i)*((mat)->max_m)+(j)]+=a)
// #define _mat_v_mul(mat,i,j,k) ((mat)->v[(i)*((mat)->max_m)+(j)]*=k)

#define _mat_swap(mat, i0, j0, i1, j1)                \
    do {                                              \
        double v = _mat_get(mat, i0, j0);             \
        _mat_set(mat, i0, j0, _mat_get(mat, i1, j1)); \
        _mat_set(mat, i1, j1, v);                     \
    } while (0)

Mat* mat_create(unsigned int m, unsigned int n) {
    Mat* mat = calloc(1, sizeof(Mat));
    mat->m = m;
    mat->max_m = m;
    mat->n = n;
    mat->max_n = n;
    mat->v = calloc(m * n, sizeof(double));
    return mat;
}

void mat_print(const Mat* mat) {
    for (int i = 0; i < mat->m; i++) {
        for (int j = 0; j < mat->n; j++)
            printf("%.2lf ", _mat_get(mat, i, j));
        printf("\n");
    }
}

void mat_input(Mat* mat) {
    for (int i = 0; i < mat->m; i++)
        for (int j = 0; j < mat->n; j++)
            scanf("%lf", mat->v + i * mat->max_m + j);
}

int mat_row_switch(Mat* mat, unsigned int i, unsigned int j) {  // return
                                                                // [error]
    if (i > mat->m || j > mat->m)
        return -1;
    if (i == j)
        return 0;
    for (int x = 0; x < mat->n; x++)
        _mat_swap(mat, i, x, j, x);
    return 0;
}

int mat_row_add(Mat* mat,
                unsigned int i,
                unsigned int j,
                double k) {  // return [error] ,R(i)+kR(j) -> R(i)
    if (i > mat->m || j > mat->m)
        return -1;
    if (k == 1.0)
        return 0;
    for (int x = 0; x < mat->n; x++)
        _mat_set(mat, i, x, _mat_get(mat, i, x) + _mat_get(mat, j, x) * k);
    return 0;
}
Mat* mat_dup(const Mat* mat) {
    Mat* dst = malloc(sizeof(Mat));
    memcpy(dst, mat, sizeof(Mat));
    double* v = malloc(mat->max_m * mat->max_n * sizeof(double));
    memcpy(v, mat->v, mat->max_m * mat->max_n * sizeof(double));
    dst->v = v;
    return dst;
}

void mat_free(Mat* mat) {
    free(mat->v);
    free(mat);
}
double mat_determinant(const Mat* _mat) {
    assert(_mat->n == _mat->m);
    Mat* mat = mat_dup(_mat);
    double k = 0.0;
    double result = 1.0;
    for (int i = 0; i < mat->n - 1; i++) {
        if (_mat_get(mat, i, i) == 0.0) {
            for (int j = i + 1; j < mat->n; j++)
                if (_mat_get(mat, j, i) != 0.0) {
                    mat_row_switch(mat, i, j);
                    result *= (j - i) % 2 == 0 ? 1 : -1;
                    break;
                }
            if (_mat_get(mat, i, i) == 0.0) {
                result = 0.0;
                goto finally;
            }
        }
        for (int j = i + 1; j < mat->n; j++) {
            k = _mat_get(mat, j, i) / _mat_get(mat, i, i) * (-1);
            mat_row_add(mat, j, i, k);
        }
        result *= _mat_get(mat, i, i);
    }
    result *= _mat_get(mat, mat->n - 1, mat->n - 1);
finally:
    mat_free(mat);
    return result;
}
double mat_cofactor(const Mat* _mat, unsigned int i0, unsigned int j0) {
    assert(_mat->m == _mat->n);
    Mat* mat = mat_create(_mat->m - 1, _mat->n - 1);
    for (int i = 0; i < mat->m; i++)
        for (int j = 0; j < mat->n; j++)
            _mat_set(
                mat, i, j,
                _mat_get(_mat, (i < i0 ? i : i + 1), (j < j0 ? j : j + 1)));
    double result =
        (i0 + j0) % 2 == 0 ? mat_determinant(mat) : -mat_determinant(mat);
    mat_free(mat);
    return result;
}
Mat* mat_adjoint(const Mat* in_mat) {
    assert(in_mat->m == in_mat->n);
    if (mat_determinant(in_mat) == 0.0)
        return NULL;
    Mat* out_mat = mat_create(in_mat->m, in_mat->n);
    for (int i = 0; i < in_mat->m; i++)
        for (int j = 0; j < in_mat->n; j++)
            _mat_set(out_mat, i, j, mat_cofactor(in_mat, i, j));
    return out_mat;
}
int main() {
    Mat* m = mat_create(3, 3);
    mat_input(m);
    mat_print(m);
    printf("\n\n");
    Mat* am = mat_adjoint(m);
    mat_print(am);
    system("pause");
}