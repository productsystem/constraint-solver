#include "raylib.h"
#include "matrix.h"
#include "gaussian_elim.h"
#include "stdlib.h"
#include "stdio.h"
#include "distance_constraint.h"
#include "math.h"
#include "constants.h"

float q[4], qdot[4], qddot[4];
float lambda[2], rhs[2];
float force[4];

void compute_qddot(float *q, float *qdot, float *qddot_out){
    float force[4] = {0};
    force[1] = g * 1.0f;
    force[3] = g * 1.0f;

    float temp_qddot[4] = {0};
    temp_qddot[1] = force[1];
    temp_qddot[3] = force[3];

    Mat *J = create_mat(2, 4);
    Mat *Jdot = create_mat(2, 4);
    Mat *W = create_mat(4, 4);
    Mat *JT = create_mat(4, 2);
    Mat *JW = create_mat(2, 4);
    Mat *JWJT = create_mat(2, 2);

    set(W, 0, 0, 1.0f); set(W, 1, 1, 1.0f); set(W, 2, 2, 1.0f); set(W, 3, 3, 1.0f);

    calc_J(J, q);
    calc_Jdot(Jdot, qdot);
    transpose_mat(J, JT);
    mult_mat(J, W, JW);
    mult_mat(JW, JT, JWJT);

    float C[2], Cdot[2], rhs[2], lambda[2];
    calc_C(C, q);
    calc_Cdot(Cdot, J, qdot);

    float *Jdotqdot = (float *)calloc(Jdot->r, sizeof(float));
    mult_mat_vec(Jdot, qdot, Jdotqdot);
    float *JWQ = (float *)calloc(JW->r, sizeof(float));
    mult_mat_vec(JW, force, JWQ);

    for (int i = 0; i < 2; i++) {
        rhs[i] = -Jdotqdot[i] - JWQ[i] - ks*C[i] - kd*Cdot[i];
    }

    if (gaussian_elim(JWJT, rhs, lambda) == -1) {
        lambda[0] = lambda[1] = 0;
    }

    float *Qcap = (float *)calloc(4, sizeof(float));
    mult_mat_vec(JT, lambda, Qcap);

    for (int i = 0; i < 4; i++) {
        qddot_out[i] = temp_qddot[i] + Qcap[i];
    }

    free(Qcap); free(Jdotqdot); free(JWQ);
    free_mat(J); free_mat(Jdot); free_mat(JT); free_mat(JW); free_mat(JWJT); free_mat(W);
}


int main(){
	InitWindow(800,600,"Constrained Dynamics");
	SetTargetFPS(120);
	for(int i = 0; i < 4; i++){
		qddot[i] = 0;
		qdot[i] = 0;
		force[i] = 0;
	}
	q[0] = l1;
	q[1] = y0;
	q[2] = l1 + l2;
	q[3] = y0;

	for(int i = 0; i < 2; i++){
		lambda[i] = 0;
		rhs[i] = 0;
	}

	Mat *J = create_mat(2,4);
	Mat *Jdot = create_mat(2,4);
	Mat *W = create_mat(4,4);
	Mat *JT = create_mat(4,2);
	Mat *JW = create_mat(2,4);
	Mat *JWJT = create_mat(2,2);

	set(W,0,0,1.0f);
	set(W,1,1,1.0f);
	set(W,2,2,1.0f);
	set(W,3,3,1.0f);
	
	while(!WindowShouldClose()){
		for (int i = 0; i < 4; i++) {
			qddot[i] = 0;
			force[i] = 0;
		}
		force[1] = g * 1.0f;
		force[3] = g * 1.0f;
		qddot[1] = force[1] * 1.0f;
		qddot[3] = force[3] * 1.0f;

		calc_J(J,q);
		calc_Jdot(Jdot,qdot);
		transpose_mat(J,JT);
		mult_mat(J,W,JW);
		mult_mat(JW,JT,JWJT);

		float *Jdotqdot = (float *)calloc(Jdot->r, sizeof(float));
		mult_mat_vec(Jdot, qdot, Jdotqdot);
		float *JWQ = (float *)calloc(JW->r, sizeof(float));
		mult_mat_vec(JW, force, JWQ);

		float C[2], Cdot[2];
		calc_C(C, q);
		calc_Cdot(Cdot, J, qdot);

		for(int i = 0; i < 2; i++){
			rhs[i] = -Jdotqdot[i] - JWQ[i] - ks*C[i] - kd*Cdot[i];
		}

		if(gaussian_elim(JWJT,rhs,lambda) == -1)
		{
			lambda[0] = 0;
			lambda[1] = 0;
		}

		float *Qcap = (float *)calloc(4, sizeof(float));
		mult_mat_vec(JT,lambda,Qcap);

		for(int i = 0; i < 4; i++){
			qddot[i] += Qcap[i];
		}

		float k1_q[4], k2_q[4], k3_q[4], k4_q[4];
		float k1_v[4], k2_v[4], k3_v[4], k4_v[4];
		float temp_q[4], temp_v[4], qddot_temp[4];

		compute_qddot(q, qdot, qddot_temp);
		for (int i = 0; i < 4; i++) {
			k1_q[i] = qdot[i];
			k1_v[i] = qddot_temp[i];
			temp_q[i] = q[i] + 0.5f * dt * k1_q[i];
			temp_v[i] = qdot[i] + 0.5f * dt * k1_v[i];
		}
		compute_qddot(temp_q, temp_v, qddot_temp);
		for (int i = 0; i < 4; i++) {
			k2_q[i] = temp_v[i];
			k2_v[i] = qddot_temp[i];
			temp_q[i] = q[i] + 0.5f * dt * k2_q[i];
			temp_v[i] = qdot[i] + 0.5f * dt * k2_v[i];
		}
		compute_qddot(temp_q, temp_v, qddot_temp);
		for (int i = 0; i < 4; i++) {
			k3_q[i] = temp_v[i];
			k3_v[i] = qddot_temp[i];
			temp_q[i] = q[i] + dt * k3_q[i];
			temp_v[i] = qdot[i] + dt * k3_v[i];
		}
		compute_qddot(temp_q, temp_v, qddot_temp);
		for (int i = 0; i < 4; i++) {
			k4_q[i] = temp_v[i];
			k4_v[i] = qddot_temp[i];
		}
		for (int i = 0; i < 4; i++) {
			q[i] += (dt / 6.0f) * (k1_q[i] + 2*k2_q[i] + 2*k3_q[i] + k4_q[i]);
			qdot[i] += (dt / 6.0f) * (k1_v[i] + 2*k2_v[i] + 2*k3_v[i] + k4_v[i]);
}

		float v1_sq = qdot[0]*qdot[0] + qdot[1]*qdot[1];
		float v2_sq = qdot[2]*qdot[2] + qdot[3]*qdot[3];
		float KE = 0.5f * (v1_sq + v2_sq);
		float PE = -g * (q[1] + q[3]);
		float E = KE + PE;

		BeginDrawing();
			ClearBackground(RAYWHITE);
			DrawCircle(x0,y0,5,BLACK);
			DrawCircle(q[0],q[1],5,BLUE);
			DrawCircle(q[2],q[3],5,RED);
			DrawLine((int)x0, (int)y0, (int)q[0], (int)q[1], DARKGRAY);
			DrawLine((int)q[0], (int)q[1], (int)q[2], (int)q[3], DARKGRAY);
			DrawText(TextFormat("Energy : %.2f", E),0,0,10,BLACK);
			float d1 = sqrtf((q[0] - x0)*(q[0] - x0) + (q[1] - y0)*(q[1] - y0));
			float d2 = sqrtf((q[2] - q[0])*(q[2] - q[0]) + (q[3] - q[1])*(q[3] - q[1]));
			DrawText(TextFormat("d1: %.2f, d2: %.2f", d1, d2), 0, 20, 10, RED);

		EndDrawing();
		free(Jdotqdot);
		free(JWQ);
		free(Qcap);
	}

	free_mat(J);
	free_mat(Jdot);
	free_mat(JT);
	free_mat(JWJT);
	free_mat(JW);
	free_mat(W);

	CloseWindow();
}