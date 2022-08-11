#include "heat.h"
#include <mpi.h>
#include <stdlib.h>

#include "communication.h"

/*
This function does following tasks:
- Each rank compute relax_jacobi and residual itself.
- After each iteration, exchange data at their boundaries with neighbors.
*/
int communicate(algoparam_t* param, MPI_Datatype column_type, MPI_Request* reqs)
{
    // Use non-blocking communication
#ifndef BLOCKING_COMMUNICATION
    return nonblocking_communication(param, column_type, reqs);
#else
    return blocking_communication(param, column_type, reqs);
#endif
}

double jacobi_iters(algoparam_t* param)
{

    // nx and ny plus ghost cells
    const int nx = param->topo_nx + 2;
    const int ny = param->topo_ny + 2;

    // Setup a datatype to send a column
    MPI_Datatype column_type;
    MPI_Type_vector(ny, 1, nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    int rank = param->rank;
    double residual = 9999999;

#ifndef BLOCKING_COMMUNICATION
    // Two request for each send, we only actually need it for non-blocking communication
    MPI_Request* reqs = (MPI_Request*) malloc(sizeof(MPI_Request) * 8);
#else
    // Else we just pass in a null pointer, it's not used anyway
    MPI_Request* reqs = NULL;
#endif

    int waitfor_count = 0;
    int finished_all = 0;

    for (int iter = 0; iter < param->maxiter; iter++) {
        // Test if we're waiting for any data
        MPI_Testall(waitfor_count, reqs, &finished_all, MPI_STATUSES_IGNORE);

        if (!finished_all) {
#ifdef COMPUTE_WHILE_WAITING
            // We didn't receive all communication so far, so do the inner part first
            double tmp_residual = inner_jacobi(param->u, param->uhelp, nx, ny, 1, 1);

            // Now wait for the rest
            MPI_Waitall(waitfor_count, reqs, MPI_STATUSES_IGNORE);

            // Now do the borders, which we left out before
            tmp_residual += outer_jacobi(param->u, param->uhelp, nx, ny, 1, 1);

            // Swap pointers now
            double* tmp = param->u;
            param->u = param->uhelp;
            param->uhelp = tmp;
#else
            // Now wait for the rest
            MPI_Waitall(waitfor_count, reqs, MPI_STATUSES_IGNORE);

            // All good all data is here just do a normal thing
            double tmp_residual = relax_jacobi(&(param->u), &(param->uhelp), nx, ny);
#endif
            // As we don't break for convergence, save the MPI_allreduce here
            residual = tmp_residual;
        } else {
            // All good all data is here just do a normal thing
            double tmp_residual = relax_jacobi(&(param->u), &(param->uhelp), nx, ny);

            // As we don't break for convergence, save the MPI_allreduce here
            residual = tmp_residual;
        }

        // Communication should happen after we're computing, we might be sending partial data in the non-blocking case
        // Returns the number of requests we should be waiting on before we can continue with our computation
        waitfor_count = communicate(param, column_type, reqs);
    }

    // Finial wait to let all communication finish
    MPI_Waitall(waitfor_count, reqs, MPI_STATUSES_IGNORE);

    // Don't forget to free the requests
#ifndef BLOCKING_COMMUNICATION
    free(reqs);
#endif
    MPI_Type_free(&column_type);
    return residual;
}

