#include "distance_constraint.h"

extern float x0;
extern float y0;
extern float l1;
extern float l2;

void calc_J(Mat *J, float *q){
    float x1 = q[0];
    float y1 = q[1];
    float x2 = q[2];
    float y2 = q[3];

    //(x1 -x0)^2 + (y1-y0)^2 - l1^2
    set(J,0,0,2*(x1-x0));
    set(J,0,1,2*(y1-y0));
    set(J,0,2,0);
    set(J,0,3,0);
    //(x2-x1)^2 + (y2-y1)^2 - l2^2
    set(J,1,0,-2*(x2-x1));
    set(J,1,1,-2*(y2-y1));
    set(J,1,2,2*(x2-x1));
    set(J,1,3,2*(y2-y1));
}

void calc_Jdot(Mat *Jdot, float *qdot){
    float dx1 = qdot[0];
    float dy1 = qdot[1];
    float dx2 = qdot[2];
    float dy2 = qdot[3];

    set(Jdot,0,0,2*dx1);
    set(Jdot,0,1,2*dy1);
    set(Jdot,0,2,0);
    set(Jdot,0,3,0);
    
    set(Jdot,1,0,-2*(dx2-dx1));
    set(Jdot,1,1,-2*(dy2-dy1));
    set(Jdot,1,2,2*(dx2-dx1));
    set(Jdot,1,3,2*(dy2-dy1));
}

void calc_C(float *C, float *q){
    float x1 = q[0], y1 = q[1];
    float x2 = q[2], y2 = q[3];
    
    C[0] = (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) - l1*l1;
    C[1] = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) - l2*l2;
}

void calc_Cdot(float *Cdot, Mat *J, float *qdot){
    mult_mat_vec(J, qdot, Cdot);
    Cdot[0] *= 0.5f;
    Cdot[1] *= 0.5f;
}