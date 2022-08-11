#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <algorithm>

#include <immintrin.h>
#include <omp.h>

#include "Stopwatch.h"

/**
 * A: M x K, ld M
 * B: K x N, ld K
 * C: M x N, ld S <<<--??
 */
void microkernel(double* A_pack, double* B_pack, double* C, int S) {

    double C_aux[M];
    int start_A, start_C;

    __assume_aligned(A_pack, ALIGNMENT);
    __assume_aligned(B_pack, ALIGNMENT);
    __assume_aligned(C, ALIGNMENT);

    for (int j = 0; j < MC/M; j++){
        for (int k = 0; k < N; k++){
            memset(C_aux, 0, M*sizeof(double));
            for (int ii = 0; ii < K; ii++){
                start_A = j*M*K + ii*M;
#pragma omp simd
                for (int jj = 0; jj < M; jj++){
                    //C_aux[j*M+k*MC+jj] += A_pack[start_A + jj] * B_pack[k*K + ii];
                    C_aux[jj] += A_pack[start_A + jj] * B_pack[k*K + ii];
                }
            }
            start_C = j*M + k*S; // <-- M  and S?
            __assume(start_C%16==0);
#pragma omp simd
            for (int jj = 0; jj < M; jj++){
                C[start_C  + jj] += C_aux[jj] ; /// <Wrong offset here!??
            }
        }
    }
}

/**
 * A: MC x K, ld S
 * B:  K x S, ld K
 * C: MC x S, ld S
 */
void GEBP(double* A, double* B_pack, double* C, double* A_pack, int S, int threadsPerTeam)
{
    // Pack A into A_pack
    // For each thread team, A_pack is initialized by its master thread
#pragma omp parallel num_threads(threadsPerTeam)
    {
#pragma omp for
        for (int j=0; j < MC/M; j++) {
            for (int ii=0; ii < K; ii++) {
                memcpy(A_pack + j*M*K + ii*M, A + j*M + ii*S , M*sizeof(double));
            }
        }
#pragma omp for nowait 
        for (int i = 0; i < S/N; i++){
            microkernel(A_pack, B_pack+i*N*K, C+i*N*S, S);
        }
    }
    /*
       double *C_aux = (double *) _mm_malloc(MC*N*sizeof(double), ALIGNMENT);
       int tmp_mc=MC,tmp_k=K, tmp_n=N;
       for (int i = 0; i < S/N; i++){
       memset(C_aux, 0, MC*N*sizeof(double));
//microkernel(A_pack, B_pack+i*N*K, C_aux, S);
try to call libxsmm_gemm() here

// Add C_aux into C
for (int j=0; j < N; j++){
#pragma omp simd
for (int k=0; k < MC; k++){
C[(i*N+j)*S + k] += C_aux[j*MC+k];
}
} // This can be re-write in the way that C_aux is MC*N, add C_aux into C can be done after all C_aux is calculated
}
_mm_free(C_aux);
     */
}

/**
 * A: S x K, ld S
 * B: K x S, ld S
 * C: S x S, ld S
 */
void GEPP(double* A, double* B, double* C, double** A_pack, double** B_pack, int S, int nTeams, int threadsPerTeam)
{

#ifdef NOFT
    // Pack B into B_pack
    for (int i=0; i < S; i++) memcpy(B_pack[0] + i*K, B + i*S , K*sizeof(double));

    int nstrides = S / MC;
    int threadTeam=omp_get_thread_num();
#pragma omp parallel for num_threads(nTeams)
    for (int i=0; i < nstrides; i++)
    {
        GEBP(A + i*MC, B_pack[0], C + i*MC, A_pack[threadTeam], S, threadsPerTeam);
    }

#else

#pragma omp parallel num_threads(nTeams)
    {
        int threadTeam=omp_get_thread_num();
        // Pack B into B_pack
#pragma omp parallel for num_threads(threadsPerTeam)
        for (int i=0; i < S; i++) memcpy(B_pack[threadTeam] + i*K, B + i*S , K*sizeof(double));

        int nstrides = S / MC;
#pragma omp for nowait
        for (int i=0; i < nstrides; i++)
        {
            GEBP(A + i*MC, B_pack[threadTeam], C + i*MC, A_pack[threadTeam], S, threadsPerTeam);
        }
    }

#endif
}

/**
 * A: S x S, ld S
 * B: S x S, ld S
 * C: S x S, ld S
 */
void GEMM(double* A, double* B, double* C, double** A_pack, double** B_pack, int S, int nTeams = 1, int threadsPerTeam = 1) {

    int nstrides = S / K;
    for (int i=0; i < nstrides; i++ )
    {
        GEPP(A + i*S*K, B + i*K, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
    }

}

