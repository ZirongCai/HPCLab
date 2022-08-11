#include <cstring>

/*  References: Fig 8 in the paper [12]
 * gebp is with assumption m_r = n_r = 1
 * gebp_re is with mr and n_r # 1
 * Assumtions:
 * - K is a multiple number of k_c
 * - M is a multiple number of m_c
 * - m_c is a multiple number of m_r
 * - N is a multiple number of n_r
 * - prefer m_r is at least equal to vector extention width (m_r>=4 for avx2)
 */

void gebp(double *subA, double *subC, double *pA, double *pB, double *C_aux)
    // Perform block*panel with Register Blocking in addition to Cache Blocking
{
    __assume_aligned(subA, ALIGNMENT);
    __assume_aligned(subC, ALIGNMENT);
    __assume_aligned(pA, ALIGNMENT);
    __assume_aligned(pB, ALIGNMENT);
    __assume_aligned(C_aux, ALIGNMENT);

    // Pack A into pA
    int nstrides_mr = m_c/m_r;
    for (int j=0; j < nstrides_mr; j++) {
        for (int ii=0; ii < k_c; ii++) {
            memcpy(pA + j*m_r*k_c + ii*m_r, subA + j*m_r + ii*M , m_r*sizeof(double));
        }
    }

    int nstrides_nr = N/n_r;
    int start_pA, start_pB, start_C;

    for (int i = 0; i < nstrides_nr; i++){
        for (int j = 0; j < nstrides_mr; j++){
            for (int k = 0; k < n_r; k++){     
                memset(C_aux, 0, m_r*sizeof(double));
                start_pB = (i*n_r + k) * k_c;

                for (int ii = 0; ii < k_c; ii++){
                    start_pA = j*m_r*k_c + ii*m_r;
#pragma omp simd
                    for (int jj = 0; jj < m_r; jj++){
                        C_aux[jj] += pA[start_pA + jj] * pB[start_pB + ii];
                    }
                }
                start_C = j*m_r + (i*n_r + k)*M; 
#pragma omp simd
                for (int jj = 0; jj < m_r; jj++){
                    subC[start_C  + jj] += C_aux[jj] ;

                }
            }
        }
    }
}
void gepp(double *subA, double *subB, double *subC, double *pA, double *pB, double *C_aux)
{
    // Pack B into pB
    for (int i=0; i < N; i++) memcpy(pB + i*k_c, subB + i*K , k_c*sizeof(double));

    int nstrides = M / m_c; 
    for (int i=0; i < nstrides; i++)
    {
        gebp(subA + i*m_c, subC + i*m_c, pA, pB, C_aux); 
    }
}

void dgemm_opt(double* A, double* B, double* C, double *pA, double *pB, double *C_aux) 
{
    int nstrides = K / k_c; 
    for (int i=0; i < nstrides; i++ )
    {
        gepp(A + i*M*k_c, B + i*k_c, C, pA, pB, C_aux); 
    }

}

