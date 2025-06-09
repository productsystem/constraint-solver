#include "matrix.h"
#include "stdlib.h"

Mat *create_mat(int r, int c){
    Mat *m = malloc(sizeof(Mat));
    m->r = r;
    m->c = c;
    m->data = calloc(r*c,sizeof(float));
    return m;
}

void free_mat(Mat *m){
    if(m){
        free(m->data);
        free(m);
    }
}

float get(Mat *m, int r, int c){
    return m->data[r * m->c + c];
}

void set(Mat *m, int r, int c, float val){
    m->data[r * m->c + c] = val;
}

void transpose_mat(Mat *m, Mat *mT){
    for (int i = 0; i < m->r; i++){
        for (int j = 0; j < m->c;j++){
            mT->data[j * mT->c + i] = m->data[i * m->c + j];
        }
    }
}

void mult_mat(Mat *a, Mat *b, Mat *out){
    for (int i = 0; i < out->r; i++) {
        for (int j = 0; j < out->c; j++) {
            float sum = 0.0f;
            for (int k = 0; k < a->c; k++) {
                sum += get(a, i, k) * get(b, k, j);
            }
            set(out, i, j, sum);
        }
    }
}

void mult_mat_vec(Mat *a, float *b, float *out){
    for (int i = 0; i < a->r; i++) {
        out[i] = 0.0f;
        for (int j = 0; j < a->c; j++) {
            out[i] += a->data[i * a->c + j] * b[j];
        }
    }
}