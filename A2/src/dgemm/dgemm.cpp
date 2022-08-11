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
#include "kernel.h"

/** Syntax: M x N, ld S
 *  A M x N sub-block of of S x something matrix
 *  in column-major layout.
 *  That is C_ij := C[j*S + i], i=0,...,M-1,  j=0,...,N-1
 */



void dgemm(double* A, double* B, double* C, int S) {
    for (int n = 0; n < S; ++n) {
        for (int m = 0; m < S; ++m) {
            for (int k = 0; k < S; ++k) {
                C[n*S + m] += A[k*S + m] * B[n*S + k];
            }
        }
    }
}

int main(int argc, char** argv) {
    int S = 4096;
    bool test = true;
    int threadsPerTeam = 1;
    int nRepeat = 10;
    if (argc <= 1) {
        printf("Usage: dgemm <S> <test> <threadsPerTeam> <repetitions>");
        return -1;
    }
    if (argc > 1) {
        S = atoi(argv[1]);
    }
    if (argc > 2) {
        test = atoi(argv[2]) != 0;
    }
    if (argc > 3) {
        threadsPerTeam = atoi(argv[3]);
    }
    if (argc > 4) {
        nRepeat = atoi(argv[4]);
    }

    omp_set_nested(1);

    int nThreads, nTeams;
#pragma omp parallel
#pragma omp master
    {
        nThreads = omp_get_num_threads(); 
    }
    threadsPerTeam = std::min(threadsPerTeam, nThreads);
    nTeams = nThreads / threadsPerTeam;

    double *A = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);
    double *B = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);
    double *C = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);
    double *A_test = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);
    double *B_test = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);
    double *C_test = (double *) _mm_malloc(S*S*sizeof(double), ALIGNMENT);

    double **A_pack = (double **) _mm_malloc(nTeams*sizeof(double), ALIGNMENT);
    double **B_pack = (double **) _mm_malloc(nTeams*sizeof(double), ALIGNMENT);
#pragma omp parallel for num_threads(nTeams)
    for (int t=0; t < nTeams; t++){
        A_pack[t] = (double *) _mm_malloc(MC*K*sizeof(double), ALIGNMENT);
        B_pack[t] = (double *) _mm_malloc(K*S*sizeof(double), ALIGNMENT);
    }


#ifdef NOFT
    for (int j = 0; j < S; ++j) {
        for (int i = 0; i < S; ++i) {
            A[j*S + i] = i + j;
            B[j*S + i] = (S-i) + (S-j);
            A[j*S + i] = 0.0;
        }
    }

#else

    int row, col;
#pragma omp parallel for num_threads(nTeams) schedule(static)
    for (int j = 0; j < S/K; ++j) {
        for (int i = 0; i < S; ++i) {
            for (int jj = 0; jj < K; ++jj) {
                col = j*K + jj;
                A[col*S + i] = i + col;
            }
        }
    }
#pragma omp parallel for num_threads(nTeams)
    for (int j = 0; j < S; ++j) {
#pragma omp parallel for num_threads(threadsPerTeam)
        for (int i = 0; i < S; ++i) {
            B[j*S + i] = (S-i) + (S-j);
        }
    }
#pragma omp parallel for num_threads(nTeams) schedule(static)
    for (int i = 0; i < S/MC; ++i) {
#pragma omp parallel for num_threads(threadsPerTeam), schedule(static) // or with chucksize=N
        for (int j = 0; j < S/N; ++j) {
            for (int ii = 0; ii < MC; ++ii) {
                for (int jj = 0; jj < N; ++jj) {
                    C[(j*N+jj)*S + (i*MC+ii)] = 0.0;
                }
            }
        }
    }
#endif
    memcpy(A_test, A, S*S*sizeof(double));
    memcpy(B_test, B, S*S*sizeof(double));
    memset(C_test, 0, S*S*sizeof(double));

    /** Check correctness of optimised dgemm */
    if (test) {
#pragma noinline
        {
            dgemm(A_test, B_test, C_test, S);
            GEMM(A, B, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
        }

        double error = 0.0;
        for (int i = 0; i < S*S; ++i) {
            double diff = C[i] - C_test[i];
            error += diff*diff;
        }
        error = sqrt(error);
        if (error > std::numeric_limits<double>::epsilon()) {
            printf("Optimised DGEMM is incorrect. Error: %e\n", error);
            //      return -1;
        }
    }

    Stopwatch stopwatch;
    double time;

    /** Test performance of microkernel */

    stopwatch.start();
    for (int i = 0; i < 10000; ++i) {
        //#pragma noinline
        microkernel(A, B, C, S);
        __asm__ __volatile__("");
    }
    time = stopwatch.stop();
    printf("Microkernel: %lf ms, %lf GFLOP/s\n", time*1.0e3, 10000*2.0*M*N*K/time * 1.0e-9);

    /** Test performance of GEBP */

    stopwatch.start();
    for (int i = 0; i < nRepeat; ++i) {
        //#pragma noinline
        GEBP(A, B, C, A_pack[0], S, threadsPerTeam);
        __asm__ __volatile__("");
    }
    time = stopwatch.stop();
    printf("GEBP: %lf ms, %lf GFLOP/s\n", time*1.0e3, nRepeat*2.0*MC*S*K/time * 1.0e-9);

    /** Test performance of optimised GEMM */

    stopwatch.start();
    for (int i = 0; i < nRepeat; ++i) {
        //#pragma noinline
        GEMM(A, B, C, A_pack, B_pack, S, nTeams, threadsPerTeam);
        __asm__ __volatile__("");
    }  
    time = stopwatch.stop();
    printf("GEMM: %lf ms, %lf GFLOP/s\n", time * 1.0e3, nRepeat*2.0*S*S*S/time * 1.0e-9);

    /** Clean up */

    _mm_free(B);
    _mm_free(C);
    _mm_free(A);
    _mm_free(A_test); _mm_free(B_test); _mm_free(C_test);
    //#pragma omp parallel num_threads(nTeams)
    for (int t = 0; t < nTeams; ++t) {
        _mm_free(A_pack[t]);
        _mm_free(B_pack[t]);
    }
    _mm_free(A_pack); _mm_free(B_pack);

    return 0;
}
