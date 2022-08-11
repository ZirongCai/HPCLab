#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <assert.h>
#include <iostream>
#include <malloc.h>

#include "Stopwatch.h"

#include "kernel.h"


void dgemm(double* A, double* B, double* C) {
    for (int n = 0; n < N; ++n) {
        for (int m = 0; m < M; ++m) {
            for (int k = 0; k < K; ++k) {
                C[n*M + m] += A[k*M + m] * B[n*K + k];
            }
        }
    }
}

int main(int argc, char** argv) {
    int repetitions = 10000;
    if (argc > 1) {
        repetitions = atoi(argv[1]);
    }

   double *A = (double *) _mm_malloc(M*K*sizeof(double), ALIGNMENT);
   double *B = (double *) _mm_malloc(K*N*sizeof(double), ALIGNMENT);
   double *C = (double *) _mm_malloc(M*N*sizeof(double), ALIGNMENT);
   double *A_test = (double *) _mm_malloc(M*K*sizeof(double), ALIGNMENT);
   double *B_test = (double *) _mm_malloc(K*N*sizeof(double), ALIGNMENT);
   double *C_test = (double *) _mm_malloc(M*N*sizeof(double), ALIGNMENT);
#ifdef OPT1
   double *pA = (double *) _mm_malloc(m_c*k_c*sizeof(double), ALIGNMENT);
   double *pB = (double *) _mm_malloc(k_c*N*sizeof(double), ALIGNMENT);
   double *C_aux = (double *) _mm_malloc(m_c*sizeof(double), ALIGNMENT);
    memset(pA, 0., m_c*k_c*sizeof(double));
    memset(pB, 0., k_c*N*sizeof(double));
#else
   double *pA, *pB, *C_aux;
#endif

   __assume_aligned(A, ALIGNMENT);
   __assume_aligned(B, ALIGNMENT);
   __assume_aligned(C, ALIGNMENT);
   __assume_aligned(A_test, ALIGNMENT);
   __assume_aligned(B_test, ALIGNMENT);
   __assume_aligned(C_test, ALIGNMENT);
   __assume_aligned(pA, ALIGNMENT);
   __assume_aligned(pB, ALIGNMENT);
   __assume_aligned(C_aux, ALIGNMENT);


    for (int j = 0; j < K; ++j) {
#pragma omp simd
        for (int i = 0; i < M; ++i) {
            A[j*M + i] = i + j;
        }
    }
    for (int j = 0; j < N; ++j) {
#pragma omp simd
        for (int i = 0; i < K; ++i) {
            B[j*K + i] = (K-i) + (N-j);
        }
    }

    memset(C, 0., M*N*sizeof(double));
    memcpy(A_test, A, M*K*sizeof(double));
    memcpy(B_test, B, K*N*sizeof(double));
    memset(C_test, 0, M*N*sizeof(double));

    /** Check correctness of optimised dgemm */
#pragma noinline
    {
        dgemm(A, B, C);
        dgemm_opt(A_test, B_test, C_test, pA, pB, C_aux);
    }
    double error = 0.0;
#pragma omp simd reduction(+: error)
    for (int i = 0; i < M*N; ++i) {
        double diff = C[i] - C_test[i];
        error += diff*diff;
    }
    error = sqrt(error);
    if (error > std::numeric_limits<double>::epsilon()) {
	    /*
        printf("Optimised DGEMM is incorrect. Error: %e\n", error);
           printf("This is A \n");
           for (int i = 0; i < M; ++i) {
           for (int j = 0; j < K; ++j) {
           std::cout << A[j*M+i] << " " ;
           }
           std::cout << std::endl; 
           }
           printf("This is B \n");
           for (int i = 0; i < K; ++i) {
           for (int j = 0; j < N; ++j) {
           std::cout << B[j*K+i] << " " ;
           }
           std::cout << std::endl; 
           }
           printf("This is the correct C \n");
           for (int i = 0; i < M; ++i) {
           for (int j = 0; j < N; ++j) {
           std::cout << C[j*M+i] << " " ;
           }
           std::cout << std::endl; 
           }
           printf("This is your C \n");
           for (int i = 0; i < M; ++i) {
           for (int j = 0; j < N; ++j) {
           std::cout << C_test[j*M+i] << " " ;
           }
           std::cout << std::endl; 
	   }
	   */
        return -1;
    }

    /** Test performance of optimised dgemm */

//#pragma noinline
 //   dgemm_opt(A, B, C);

    Stopwatch stopwatch;
    stopwatch.start();

#pragma noinline
    for (int r = 0; r < repetitions; ++r) {
   //     dgemm(A, B, C);
        dgemm_opt(A_test, B_test, C_test, pA, pB, C_aux);
        //__asm__ __volatile__("");
    }

    double time = stopwatch.stop();

#ifndef PRINT_FORMAT
    printf("%lf ms, %lf GFLOP/s\n", time * 1.0e3, repetitions*2.0*M*N*K/time * 1.0e-9);
#else
#ifdef OPT1
    printf("K=%d   M=%d   N=%d   m_c=%d   k_c=%d   m_r=%d   n_r=%d \t %lf GFLOPs\n", K, M, N, m_c, k_c, m_r, n_r, repetitions*2.0*M*N*K/time * 1.0e-9);
#else
    printf("K=%d   M=%d   N=%d   block_size=%d \t %lf GFLOPs\n", K, M, N, BLOCK_SIZE, repetitions*2.0*M*N*K/time * 1.0e-9);
#endif
#endif

    /** Clean up */
/*
    free(A); free(B); free(C);
    free(A_test); free(B_test); free(C_test);
    free(pA); free(pB); free(C_aux);
*/
    _mm_free(A); _mm_free(B); _mm_free(C);
    _mm_free(A_test); _mm_free(B_test); _mm_free(C_test);
#ifdef OPT1
    _mm_free(pA); _mm_free(pB); _mm_free(C_aux);
#endif

    return 0;
}
