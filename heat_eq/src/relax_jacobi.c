/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"

#define COMP(a, b) (0.25 * (u[(a) + ((b) -1)] + u[(a) + ((b) + 1)] + u[(iim1) + (b)] + u[(iip1) + (b)]))
#define SUM()                           \
    {                                   \
        double diff = unew - u[ii + j]; \
        sum += diff * diff;             \
    }

double outer_jacobi(double* u, double* utmp, unsigned sizex, unsigned sizey, unsigned borderx, unsigned bordery)
{
    double sum = 0.0;

    for (int i = 1; i < sizey - 1; ++i) {
        int ii = i * sizex;
        int iim1 = (i - 1) * sizex;
        int iip1 = (i + 1) * sizex;

        if (i < 1 + bordery || i > sizey - 2 - bordery) {
            // Process entire row
            for (int j = 1; j < sizex - 1; ++j) {

                double unew = COMP(ii, j);
                SUM();
                #pragma vector nontemporal
                utmp[ii + j] = unew;
            }
        } else {
            // Process first and last elements
            for (int jj = 0; jj < bordery; ++jj) {
                int j = 1 + jj;

                double unew = COMP(ii, j);
                SUM();
                #pragma vector nontemporal
                utmp[ii + j] = unew;
            }

            for (int jj = 0; jj < bordery; ++jj) {
                int j = sizex - 2 - jj;

                double unew = COMP(ii, j);
                SUM();
                #pragma vector nontemporal
                utmp[ii + j] = unew;
            }
        }
    }

    return (sum);
}

/// Normal jacoby, but don't work on [1, sizex[ and [1, sizey[, but rather in [1 + offsetx, sixe - offsetx[
/// and [1 + offsety, sizey - offsety[. Basically this shrinks the area by offsetx and offsety
double inner_jacobi(double* u, double* utmp, unsigned sizex, unsigned sizey, unsigned offsetx, unsigned offsety)
{
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum) schedule(guided)
    for (int i = 1 + offsety; i < sizey - 1 - offsety; i++) {
        int ii = i * sizex;
        int iim1 = (i - 1) * sizex;
        int iip1 = (i + 1) * sizex;
#pragma omp simd reduction(+ : sum)
        for (int j = 1 + offsetx; j < sizex - 1 - offsety; j++) {
            double unew = COMP(ii, j);

            SUM();

            #pragma vector nontemporal
            utmp[ii + j] = unew;
        }
    }

    return (sum);
}

double relax_jacobi(double** u1, double** utmp1, unsigned sizex, unsigned sizey)
{
    double* utmp = *utmp1;
    double* u = *u1;

    double sum = inner_jacobi(u, utmp, sizex, sizey, 0, 0);

    // Do a swap
    *u1 = utmp;
    *utmp1 = u;

    return (sum);
}
