#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <stdio.h>
#include <cstring>

#include "header.h"

void timer_start(struct timeval &begin)
{
    gettimeofday(&begin,(struct timezone *)0);
}

double timer_stop(struct timeval &begin, struct timeval &end)
{
    gettimeofday(&end,(struct timezone *)0);
    double seconds, useconds;
    double ret, tmp;

    if (end.tv_usec >= begin.tv_usec)
    {
        seconds = (double)end.tv_sec - (double)begin.tv_sec;
        useconds = (double)end.tv_usec - (double)begin.tv_usec;
    }
    else
    {
        seconds = (double)end.tv_sec - (double)begin.tv_sec;
        seconds -= 1;                   // Correction
        useconds = (double)end.tv_usec - (double)begin.tv_usec;
        useconds += 1000000;            // Correction
    }

    // get time in seconds
    tmp = (double)useconds;
    ret = (double)seconds;
    tmp /= 1000000;
    ret += tmp;

    return ret;
}
/**
 * stores a given grid into a file
 *
 * @param grid the grid that should be stored
 * @param filename the filename
 */
void store_grid(struct_param *param, double* grid, std::string filename)
{
    int start_offset_x = 0;
    int end_offset_x = 0;
    int start_offset_y = 0;
    int end_offset_y = 0;

    auto mesh_width = param->mesh_width;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;


    // if the rank contains boundary cells, incl boundary into the iteration range
    if (param->topo_bottom == MPI_PROC_NULL) { // No bottom neighbor
        start_offset_y = -1;
    }
    if (param->topo_top == MPI_PROC_NULL) { // No top neighbor
        end_offset_y = 1;
    }
    if (param->topo_left == MPI_PROC_NULL) { // No left neighbor
        start_offset_x = -1;
    }
    if (param->topo_right == MPI_PROC_NULL) { // No right neighbor
        end_offset_x = 1;   
    }

    int global_i, global_j;
    global_j = param->topo_x * param->avrg_nx;
    global_i = param->topo_y * param->avrg_ny;

    // Each rank write to file one in order
    for (int coordy = 0; coordy < param->dimy; coordy++){
        for (int coordx = 0; coordx < param->dimx; coordx++){
            if ((param->topo_x == coordx) && (param->topo_y == coordy)){
                std::fstream filestr;
                filestr.open (filename.c_str(), std::fstream::out | std::fstream::app);

                for (int i = 1 + start_offset_y; i < ny - 1 + end_offset_y; i++){
                    for (int j = 1 + start_offset_x; j < nx - 1 + end_offset_x; j++)
                    {
                        filestr << mesh_width*(global_i + i) << " " << mesh_width*(global_j + j) << " " << grid[(i*nx+j)] << std::endl;
                    }
                    //filestr << std::endl;
                }
                filestr.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}
/**
 * calculate the grid's initial values for given grid points
 *
 * @param x the x-coordinate of a given grid point
 * @param y the y-coordinate of a given grid point
 *
 * @return the initial value at position (x,y)
 */
double eval_init_func(double x, double y)
{
    return (x*x)*(y*y);
}

/**
 * initializes a given grid: inner points are set to zero
 * boundary points are initialized by calling eval_init_func
 *
 * @param grid the grid to be initialized
 */
void init_grid(struct_param *param)
{
    // set all points to zero
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;
    auto mesh_width = param->mesh_width;
    /*
       for (int i = 0; i < nx*ny ; i++)
       {
       grid[i] = 0.0;
       }
     */
    memset(param->grid, 0.0, nx*ny*sizeof(double));

    int global_i, global_j;

    // x - boundaries
    /* Remove since eval_init_func return 0.0
       if (param->topo_bottom == MPI_PROC_NULL) { // No bottom neighbor
       for (int i = 0; i < nx; i++)
       grid[i] = eval_init_func(0.0, ((double)i)*mesh_width);
       }
     */
    if (param->topo_top == MPI_PROC_NULL) { // No top neighbor
        global_i = param->topo_x * param->avrg_nx;
        for (int i = 0; i < nx; i++)
            param->grid[i + nx*(ny-1)] = eval_init_func(1.0, ((double)global_i + i)*mesh_width);
    }

    // y-boundaries
    /* Remove since eval_init_func return 0.0
       if (param->topo_left == MPI_PROC_NULL) { // No left neighbor
       for (int i = 0; i < ny; i++)
       grid[i*nx] = eval_init_func(((double)i)*mesh_width, 0.0);
       }
     */
    if (param->topo_right == MPI_PROC_NULL) { // No right neighbor
        global_j = param->topo_y * param->avrg_ny;
        for (int j = 1; j < ny-1; j++)
            param->grid[(j*nx) + (nx-1)] = eval_init_func(((double) global_j + j)*mesh_width, 1.0);
    }
}

/**
 * initializes the right hand side, we want to keep it simple and
 * solve the Laplace equation instead of Poisson (-> b=0)
 *
 * @param b the right hand side
 */
void init_b(struct_param* param)
{
    // set all points to zero
    /*
       for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
       {
       b[i] = 0.0;
       }
     */
    memset(param->b, 0.0, param->topo_nx * param->topo_ny * sizeof(double));
}

/**
 * copies data from one grid to another
 *
 * @param dest destination grid
 * @param src source grid
 */
void g_copy(double* dest, double* src, std::size_t grid_points_1d)
{
    for (int i = 0; i < grid_points_1d*grid_points_1d; i++)
    {
        dest[i] = src[i];
    }
}

/**
 * calculates the dot product of the two grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid1 first grid
 * @param grid2 second grid
 */
double g_dot_product(struct_param* param, double* grid1, double* grid2)
{
    double tmp = 0.0;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    for (int i = 1; i < ny-1; i++)
    {
        for (int j = 1; j < nx-1; j++)
        {
            tmp += (grid1[(i*nx)+j] * grid2[(i*nx)+j]);
        }
    }

    return tmp;
}

/**
 * scales a grid by a given scalar (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param grid grid to be scaled
 * @param scalar scalar which is used to scale to grid
 */
void g_scale(struct_param *param, double* grid, double scalar)
{
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    for (int i = 1; i < ny-1; i++)
    {
        for (int j = 1; j < nx-1; j++)
        {
            grid[(i*nx)+j] *= scalar;
        }
    }
}

/**
 * implements BLAS's Xaxpy operation for grids (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 *
 * @param dest destination grid
 * @param src source grid
 * @param scalar scalar to scale to source grid
 */
void g_scale_add(struct_param *param, double* dest, double* src, double scalar)
{
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    for (int i = 1; i < ny-1; i++)
    {
        for (int j = 1; j < nx-1; j++)
        {
            dest[(i*nx)+j] += (scalar*src[(i*nx)+j]);
        }
    }
}

/**
 * implements the the 5-point finite differences stencil (only inner grid points are modified due 
 * to Dirichlet boundary conditions)
 * 
 * @param grid grid for which the stencil should be evaluated
 * @param result grid where the stencil's evaluation should be stored
 */
void g_product_operator(struct_param* param, double* grid, double* result)
{
    auto mesh_width = param->mesh_width;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    for (int i = 1; i < ny - 1; i++)
    {
        for (int j = 1; j < nx - 1; j++)
        {
            result[i * nx + j] =  (
                    (4.0 * grid[i * nx + j]) 
                    - grid[(i+1) * nx + j]
                    - grid[(i-1) * nx + j]
                    - grid[ i    * nx + j+1]
                    - grid[ i    * nx + j-1]
                    ) * (mesh_width*mesh_width);
        }
    }
}

