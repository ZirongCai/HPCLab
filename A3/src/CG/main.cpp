/******************************************************************************
 * Copyright (C) 2011 Technische Universitaet Muenchen                         *
 * This file is part of the training material of the master's course           *
 * Scientific Computing                                                        *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <x86intrin.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "header.h"

/**
 * main application
 *
 * @param argc number of cli arguments
 * @param argv values of cli arguments
 */
/// store begin timestep
struct timeval begin;
/// store end timestep
struct timeval end;


int main(int argc, char* argv[])
{
    //global variable!

    // param stores all data the rank needs
    struct_param param;

    // check if all parameters are specified
    if (argc != 6)
    {
        std::cout << std::endl;
        std::cout << "meshwidth" << std::endl;
        std::cout << "cg_max_iterations" << std::endl;
        std::cout << "cg_eps" << std::endl;
        std::cout << "mpi_topo_x" << std::endl;
        std::cout << "mpi_topo_y" << std::endl;
        std::cout << std::endl;
        std::cout << "example:" << std::endl;
        std::cout << "./app 0.125 100 0.0001 1 1" << std::endl;
        std::cout << std::endl;

        return -1;
    }

    // read cli arguments
    double mesh_width = atof(argv[1]);
    size_t cg_max_iterations = atoi(argv[2]);
    double cg_eps = atof(argv[3]);
    int dim[2];
    dim[0] = atoi(argv[4]);
    dim[1] = atoi(argv[5]);

    // --------------------------------------------------------------------------------
    // Initialize MPI and param
    // --------------------------------------------------------------------------------
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //int period[2], reorder;
    //int coord[2], id;

    // Check if provided topology matches with mpi_size
    if (size != (dim[0] * dim[1])){
        fprintf(stderr, "Can not work with the provided topology.\n");
        fprintf(stderr, "Topo_x by topo_y should be the number of mpi processes.\n");
        MPI_Finalize();
        return 0;
    }

    // Set rank
    param.rank = rank;
    param.mpi_size = size;
    param.dimx = dim[0];
    param.dimy = dim[1];

    int period[2];
    period[0] = 0;
    period[1] = 0;

    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &param.comm_2d);
    if (rank == 0) {
        int st = 0;
        MPI_Topo_test(param.comm_2d, &st);
        if (st == MPI_CART) {
            fprintf(stderr, "This test uses 2D Cartesian topology! \n");
        }
    }

    // Finding coordinators of each rank
    int coord[2];
    MPI_Cart_coords(param.comm_2d, rank, 2, coord);

    // Finding the neighbors
    int right = -1;
    int left = -1;
    int top = -1;
    int bottom = -1;
    MPI_Cart_shift(param.comm_2d, 0, 1, &left, &right);
    MPI_Cart_shift(param.comm_2d, 1, -1, &top, &bottom);

    // Write data in to param struct
    param.topo_x = coord[0];
    param.topo_y = coord[1];

    param.topo_bottom = bottom;
    param.topo_top = top;
    param.topo_left = left;
    param.topo_right = right;

    //------------------------------------------------------------------------------------
    // Data initialization
    //------------------------------------------------------------------------------------
    // calculate grid points per dimension
    param.mesh_width = mesh_width;
    std::size_t	grid_points_1d = (std::size_t)(1.0/mesh_width)+1;
    auto actual_gridsize = grid_points_1d - 2;

    param.grid_points_1d = grid_points_1d;
    param.actual_gridsize = actual_gridsize;

    // Calculate local sizes, don't consider boundary cells yet
    std::size_t tmp_nx = ((actual_gridsize + param.dimx - 1) / param.dimx);
    std::size_t tmp_ny = ((actual_gridsize + param.dimy - 1) / param.dimy);
    param.avrg_nx = tmp_nx;
    param.avrg_ny = tmp_ny;

    // find start and end coord in domain for each block in topology
    // start at 1, add nx * topo_x then subtract 1 for ghost cell at the left/top
    std::size_t startx = 1 + (param.topo_x * tmp_nx) - 1;
    std::size_t starty = 1 + (param.topo_y * tmp_ny) - 1;

    // Add two at the end, as wee start at 1, plus 1 ghost cell
    std::size_t endx = std::min((param.topo_x + 1) * tmp_nx + 2, grid_points_1d);
    std::size_t endy = std::min((param.topo_y + 1) * tmp_ny + 2, grid_points_1d);

    param.glob_startx = startx;
    param.glob_starty = starty;

    param.glob_endx = endx;
    param.glob_endy = endy;

    // store the sizes which including ghosh layers
    param.topo_nx = endx - startx;
    param.topo_ny = endy - starty ;

    // allocate the grid and rights hand side locally
    param.grid = (double*) _mm_malloc( param.topo_nx * param.topo_ny * sizeof(double), 64);
    param.b 	 = (double*) _mm_malloc( param.topo_nx * param.topo_ny * sizeof(double), 64);

    init_grid(&param);
    store_grid(&param , param.grid, "initial_condition_parallel.gnuplot");
    init_b(&param);
    store_grid(&param, param.b,  "b_parallel.gnuplot");

    // solve Poisson equation using CG method
    timer_start(begin);
    solve(&param, cg_max_iterations, cg_eps);
    double time = timer_stop(begin, end);
    store_grid(&param, param.grid, "solution_parallel.gnuplot");

    if (rank==0) std::cout << std::endl << "Needed time: " << time << " s" << std::endl << std::endl;

    _mm_free(param.grid);
    _mm_free(param.b);

    MPI_Finalize();

    return 0;
}
