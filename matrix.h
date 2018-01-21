#ifndef _matrix_h
#define _matrix_h

typedef struct Mat {
    unsigned int m, n;
    unsigned int max_m, max_n;
    double* v;
} Mat;

Mat* mat_create(unsigned int m, unsigned int n);
void mat_print(const Mat* mat);
void mat_input(Mat* mat);
int mat_row_switch(Mat* mat, unsigned int i, unsigned int j);
int mat_row_add(Mat* mat, unsigned int i, unsigned int j, double k);
#define mat_row_mul(mat, r, k) (mat_row_add(mat, r, k))
Mat* mat_dup(const Mat* mat);
void mat_free(Mat* mat);
double mat_determinant(const Mat* _mat);
double mat_cofactor(const Mat* _mat, unsigned int i0, unsigned int j0);
Mat* mat_adjoint(const Mat* in_mat);
#endif