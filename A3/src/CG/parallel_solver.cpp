#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include <cstring>
#include <mpi.h>
#include "header.h"

void exchange_ghoshlayers(struct_param* param, double *grid){
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    // Create some tags for communication
    int toptag = 0;
    int downtag = 1;
    int lefttag = 2;
    int righttag = 3;

    // Exchange top most row to the top neighbors:
    MPI_Sendrecv(grid + (ny-2)*nx, nx, MPI_DOUBLE, param->topo_top, toptag, 
            grid, nx, MPI_DOUBLE, param->topo_bottom, toptag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange row most row to the bottom neighbors:
    MPI_Sendrecv(grid + nx, nx, MPI_DOUBLE, param->topo_bottom, downtag,
            grid + (ny-1)*nx, nx, MPI_DOUBLE, param->topo_top, downtag, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Data type for sending columns
    MPI_Datatype column_type;
    MPI_Type_vector(ny, 1, nx, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    // Exchange left most column to the left neighbors:
    MPI_Sendrecv(grid + 1, 1, column_type, param->topo_left, lefttag, grid + nx - 1, 1, column_type,
            param->topo_right, lefttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange right most column to the right neighbors:
    MPI_Sendrecv(grid + nx - 2, 1, column_type, param->topo_right, righttag, grid , 1, column_type,
            param->topo_left, righttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

}

void solve(struct_param *param, std::size_t cg_max_iterations, double cg_eps)
{
    if (param->rank == 0) std::cout << "Starting Conjugated Gradients" << std::endl;

    double eps_squared = cg_eps*cg_eps;
    std::size_t needed_iters = 0;

    auto nx = param->topo_nx;
    auto ny = param->topo_ny;
    auto grid = param->grid;
    auto b = param -> b;

    // define temporal vectors
    double* q = (double*)_mm_malloc(nx*ny*sizeof(double), 64);
    double* r = (double*)_mm_malloc(nx*ny*sizeof(double), 64);
    double* d = (double*)_mm_malloc(nx*ny*sizeof(double), 64);
    double* b_save = (double*)_mm_malloc(nx*ny*sizeof(double), 64);

    memcpy(q, grid, nx*ny*sizeof(double));
    memcpy(r, grid, nx*ny*sizeof(double));
    memcpy(d, grid, nx*ny*sizeof(double));
    memcpy(b_save, b, nx*ny*sizeof(double));

    double delta_0 = 0.0;
    double delta_old = 0.0;
    double delta_new = 0.0;
    double beta = 0.0;
    double a = 0.0;
    double residuum = 0.0;

    g_product_operator(param, grid, d);
    g_scale_add(param, b, d, -1.0);
    memcpy(r, b, nx*ny*sizeof(double));
    memcpy(d, r, nx*ny*sizeof(double));

    // calculate starting norm
    auto delta_local = g_dot_product(param, r, r);
    double delta_global_new, delta_global_old;

    MPI_Allreduce(& delta_local, &delta_global_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delta_0 = delta_global_new*eps_squared;

    residuum = (delta_0/eps_squared);

    if (param->rank == 0){
        std::cout << "Starting norm of residuum: " << (delta_0/eps_squared) << std::endl;
        std::cout << "Target norm:               " << (delta_0) << std::endl;
    }

    while ((needed_iters < cg_max_iterations) && (delta_global_new > delta_0))
    {
        // q = A*d
        exchange_ghoshlayers(param,d);
        g_product_operator(param, d, q);

        // a = d_new / d.q
        auto a_local = g_dot_product(param, d, q);
        MPI_Allreduce(&a_local, &a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        a = delta_global_new/a;

        // x = x + a*d
        g_scale_add(param, grid, d, a);

        if ((needed_iters % 50) == 0)
        {
            memcpy(b, b_save , nx*ny*sizeof(double));
            exchange_ghoshlayers(param,grid);
            g_product_operator(param, grid, q);
            g_scale_add(param, b, q, -1.0);
            memcpy(r, b , nx*ny*sizeof(double));
        }
        else
        {
            // r = r - a*q
            g_scale_add(param, r, q, -a);
        }

        // calculate new deltas and determine beta
        delta_global_old = delta_global_new;
        delta_local = g_dot_product(param, r, r);
        MPI_Allreduce(&delta_local, &delta_global_new, 1, MPI_DOUBLE, MPI_SUM, param->comm_2d);

        beta = delta_global_new/delta_global_old;

        // adjust d
        g_scale(param, d, beta);
        g_scale_add(param, d, r, 1.0);

        residuum = delta_global_new;
        needed_iters++;

        if (param->rank == 0){
            std::cout << "(iter: " << needed_iters << ")delta: " << delta_global_new << std::endl;
        }

    }

    if (param->rank == 0){
        std::cout << "Number of iterations: " << needed_iters << " (max. " << cg_max_iterations << ")" << std::endl;
        std::cout << "Final norm of residuum: " << delta_global_new << std::endl;
    }

    _mm_free(d);
    _mm_free(q);
    _mm_free(r);
    _mm_free(b_save);
}
