#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <stdio.h>
#include <cstring>
#include <immintrin.h>

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

	std::stringstream buffer;
	MPI_File file;
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);

	for (int i = 1 + start_offset_y; i < ny - 1 + end_offset_y; i++){
		for (int j = 1 + start_offset_x; j < nx - 1 + end_offset_x; j++)
		{
			buffer << mesh_width*(global_i + i) << " " << mesh_width*(global_j + j) << " " << grid[i*nx+j] << std::endl;
			//filestr << mesh_width*(global_i + i) << " " << mesh_width*(global_j + j) << " " << grid[(i*nx+j)] << std::endl;
		}
		buffer << std::endl;
	}
	MPI_File_write_ordered(file, buffer.str().c_str(), buffer.str().length(), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&file);
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
    auto grid = param->grid;
    auto mesh_width = param->mesh_width;

    __assume_aligned(grid, 64);

#pragma omp parallel for
    for (int i = 0; i < ny; i++)
    {
        for (int j = 1; j < nx; j++)
        {
            grid[i*nx+j] = 0.0;
        }
    }

    int global_i, global_j;

    if (param->topo_top == MPI_PROC_NULL) { // No top neighbor
        global_i = param->topo_x * param->avrg_nx;
        for (int i = 0; i < nx; i++)
            grid[i + nx*(ny-1)] = eval_init_func(1.0, ((double)global_i + i)*mesh_width);
    }

    if (param->topo_right == MPI_PROC_NULL) { // No right neighbor
        global_j = param->topo_y * param->avrg_ny;
        for (int j = 1; j < ny-1; j++)
            grid[(j*nx) + (nx-1)] = eval_init_func(((double) global_j + j)*mesh_width, 1.0);
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
    auto b = param->b;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    __assume_aligned(b, 64);
#pragma omp parallel for
    for (int i = 0; i < ny; i++)
    {
        for (int j = 1; j < nx; j++)
        {
            b[i*nx+j] = 0.0;
        }
    }
}
void g_copy(struct_param* param, double* dest, double* src)
{
    __assume_aligned(dest, 64);
    __assume_aligned(src, 64);
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;
#pragma omp parallel for
    for (int i = 0; i < ny; i++)
    {
        for (int j = 1; j < nx; j++)
        {
            dest[i*nx+j] = src[i*nx+j] ;
        }
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
    __assume_aligned(grid1, 64);
    __assume_aligned(grid2, 64);

    double tmp = 0.0;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;
#pragma omp parallel for reduction(+:tmp)
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
    __assume_aligned(grid, 64);
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

#pragma omp parallel for
    for (int i = 1; i < ny-1; i++)
    {
#pragma ivdep
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
    __assume_aligned(dest, 64);
    __assume_aligned(src, 64);
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

#pragma omp parallel for
    for (int i = 1; i < ny-1; i++)
    {
#pragma ivdep
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
    __assume_aligned(grid, 64);
    __assume_aligned(result, 64);


    auto mesh_width = param->mesh_width;
    auto nx = param->topo_nx;
    auto ny = param->topo_ny;

    __assume(nx%16==0);

#pragma omp parallel for
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

