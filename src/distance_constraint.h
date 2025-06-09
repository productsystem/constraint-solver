#pragma once
#include "matrix.h"

void calc_J(Mat *J, float *q);
void calc_Jdot(Mat *Jdot, float *qdot);
void calc_C(float *C, float *q);
void calc_Cdot(float *Cdot, Mat *J, float *qdot);
