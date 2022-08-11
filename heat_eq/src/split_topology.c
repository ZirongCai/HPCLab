#include "heat.h"
#include <mpi.h>
/* This function does following tasks:
   - Take the input grid and create virtual topology
   - Stored coordinator and neighbors for each rank
 */

void split_topography(algoparam_t* param, MPI_Comm* comm)
{
    // Set the topography at a default
    param->topo_x = 0;
    param->topo_y = 0;

    param->topo_down = -1;
    param->topo_top = -1;
    param->topo_left = -1;
    param->topo_right = -1;

    int rank = param->rank;
    int act_size = param->dimx * param->dimy;

    if (act_size == 1) {
        // With only a sizeof one, we just run single threaded (basically no MPI)
        fprintf(stderr, "The following test is executed with a single process (without MPI) \n");

        // Write data in to algo struct
        param->topo_x = 0;
        param->topo_y = 0;

        param->topo_down = 0;
        param->topo_top = 0;
        param->topo_left = 0;
        param->topo_right = 0;
    } else {
        // Passing value for dims
        int dim[2];
        dim[0] = param->dimx;
        dim[1] = param->dimy;

        int period[2];
        period[0] = 0;
        period[1] = 0;

        int reorder = 1;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, comm);

        if (rank == ROOT_RANK) {
            int st = 0;
            MPI_Topo_test(*comm, &st);
            if (st == MPI_CART) {
                fprintf(stderr, "This test uses Cartesian topology! \n");
            }
        }

        MPI_Barrier(*comm);

        // Finding coordinators of each rank
        int coord[2];
        MPI_Cart_coords(*comm, rank, 2, coord);

        // Finding the neighbors
        int right = -1;
        int left = -1;
        int top = -1;
        int bottom = -1;
        MPI_Cart_shift(*comm, 0, 1, &left, &right);
        MPI_Cart_shift(*comm, 1, 1, &top, &bottom);

        // Write data in to algo struct
        param->topo_x = coord[0];
        param->topo_y = coord[1];

        param->topo_down = bottom;
        param->topo_top = top;
        param->topo_left = left;
        param->topo_right = right;
    }
}
