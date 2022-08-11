#include <algorithm>
void dgemm_opt(double* A, double* B, double* C, double* nouse1, double* nouse2, double* nouse3) {

    __assume_aligned(A, ALIGNMENT); 
    __assume_aligned(B, ALIGNMENT); 
    __assume_aligned(C, ALIGNMENT);   

    int ii, kk, jj, m, n, k, B_reused;

    for(n = 0; n < N; n += BLOCK_SIZE){
        for(k = 0; k < K; k += BLOCK_SIZE){
            for(m = 0; m < M; m += BLOCK_SIZE){
                for(jj = n; jj < n+BLOCK_SIZE; jj++){
                    for(kk = k; kk < k+BLOCK_SIZE; kk++){
                        B_reused = B[jj*K+kk];
                        for(ii = m; ii < m+BLOCK_SIZE; ii++){
                            //A_reused = A[K*ii+jj];
                            C[jj*M+ii] += A[kk*M+ii]* B_reused;
                        }
                    }
                }
            }
        }
    }
}
/*
void dgemm_opt(double* A, double* B, double* C, double* nouse1, double* nouse2, double* nouse3) {

    __assume_aligned(A, ALIGNMENT); 
    __assume_aligned(B, ALIGNMENT); 
    __assume_aligned(C, ALIGNMENT);   

    int ii, kk, jj, m, n, k, B_reused;

    for(n = 0; n < N; n += BLOCK_SIZE){
        for(k = 0; k < K; k += BLOCK_SIZE){
            for(m = 0; m < M; m += BLOCK_SIZE){
                for(jj = n; jj < n+BLOCK_SIZE; jj++){
                    for(kk = k; kk < k+BLOCK_SIZE; kk++){
                        B_reused = B[jj*K+kk];
                        for(ii = m; ii < m+BLOCK_SIZE; ii+=16){
                            //A_reused = A[K*ii+jj];
                            C[jj*M+ii] += A[kk*M+ii]* B_reused;
                            C[jj*M+ii+1] += A[kk*M+ii+1]* B_reused;
                            C[jj*M+ii+2] += A[kk*M+ii+2]* B_reused;
                            C[jj*M+ii+3] += A[kk*M+ii+3]* B_reused;
                            C[jj*M+ii+4] += A[kk*M+ii+4]* B_reused;
                            C[jj*M+ii+5] += A[kk*M+ii+5]* B_reused;
                            C[jj*M+ii+6] += A[kk*M+ii+6]* B_reused;
                            C[jj*M+ii+7] += A[kk*M+ii+7]* B_reused;
                            C[jj*M+ii+8] += A[kk*M+ii+8]* B_reused;
                            C[jj*M+ii+9] += A[kk*M+ii+9]* B_reused;
                            C[jj*M+ii+10] += A[kk*M+ii+10]* B_reused;
                            C[jj*M+ii+11] += A[kk*M+ii+11]* B_reused;
                            C[jj*M+ii+12] += A[kk*M+ii+12]* B_reused;
                            C[jj*M+ii+13] += A[kk*M+ii+13]* B_reused;
                            C[jj*M+ii+14] += A[kk*M+ii+14]* B_reused;
                            C[jj*M+ii+15] += A[kk*M+ii+15]* B_reused;
                        }
                    }
                }
            }
        }
    }
}
*/
