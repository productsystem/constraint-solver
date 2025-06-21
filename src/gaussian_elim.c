#include "gaussian_elim.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

//NOT GAUSSIAN ELIM IT IS CG METHOD

float dot(const float *a, const float *b, int n) {
    float s = 0.0f;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

void mat_vec_mul(const Mat *A, const float *x, float *out) {
    int n = A->r;
    for (int i = 0; i < n; i++) {
        out[i] = 0.0f;
        for (int j = 0; j < n; j++) {
            out[i] += get(A, i, j) * x[j];
        }
    }
}

void vec_add(float *out, const float *a, const float *b, float alpha, int n) {
    for (int i = 0; i < n; i++) {
        out[i] = a[i] + alpha * b[i];
    }
}

void vec_sub(float *out, const float *a, const float *b, int n) {
    for (int i = 0; i < n; i++) {
        out[i] = a[i] - b[i];
    }
}

int conjugate_gradient(const Mat *A, const float *b, float *x, float tol, int max_iter) {
    int n = A->r;
    float *r = malloc(n * sizeof(float));
    float *p = malloc(n * sizeof(float));
    float *Ap = malloc(n * sizeof(float));
    if (!r || !p || !Ap) return -1;

    mat_vec_mul(A, x, r);
    for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
    memcpy(p, r, n * sizeof(float));

    float rsold = dot(r, r, n);

    for (int i = 0; i < max_iter; i++) {
        mat_vec_mul(A, p, Ap);
        float denom = dot(p, Ap, n);
        if (fabsf(denom) < 1e-8f) return -1;
        float alpha = rsold / denom;

        for (int j = 0; j < n; j++) x[j] += alpha * p[j];

        for (int j = 0; j < n; j++) r[j] -= alpha * Ap[j];

        float rsnew = dot(r, r, n);
        if (sqrtf(rsnew) < tol) break;

        for (int j = 0; j < n; j++) p[j] = r[j] + (rsnew / rsold) * p[j];
        rsold = rsnew;
    }

    free(r);
    free(p);
    free(Ap);
    return 0;
}

int solve(Mat *A, float *b, float *out) {
    memset(out, 0, sizeof(float) * A->r); 
    return conjugate_gradient(A, b, out, 1e-5f, 1000);
}

