#pragma once

typedef struct{
    int r,c;
    float *data;
} Mat;

Mat *create_mat(int r, int c);
void free_mat(Mat *m);
float get(Mat *m, int r, int c);
void set(Mat *m, int r, int c, float val);
void transpose_mat(Mat *m, Mat *mT);
void mult_mat(Mat *a, Mat *b, Mat *out);
void mult_mat_vec(Mat *a, float *b, float *out);