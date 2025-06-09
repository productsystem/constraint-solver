#include "gaussian_elim.h"
#include "stdlib.h"
#include "math.h"

int gaussian_elim(Mat *A, float *b, float *out){
    int n = A->r;
    for(int k = 0; k <n; k++){
        int max_r = k;
        for(int i= k+1;i<n;i++){
            if(fabs(get(A,i,k)) > fabs(get(A,max_r,k))){
                max_r = i;
            }
        }

        for(int j = 0; j < n;j++){
            float t = get(A, k,j);
            set(A,k,j,get(A,max_r,j));
            set(A,max_r,j,t);
        }

        float temp = b[k];
        b[k] = b[max_r];
        b[max_r] = temp;
        if (get(A, k, k) == 0.0f){
            return -1;
        }

        for(int i = k + 1; i< n;i++){
            float mult = get(A,i,k) / get(A,k,k);
            for(int j = k; j < n; j++){
                set(A,i,j,(get(A,i,j) - mult * get(A,k,j)));
            }
            b[i] -= mult * b[k];
        }
    }

    for(int i = n-1; i>=0;i--){
        float s = b[i];
        for(int j = i+1;j<n;j++){
            s-=get(A,i,j) * out[j];
        }
        out[i] = s/get(A,i,i);
    }

    return 0;
}