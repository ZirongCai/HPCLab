
#include "communication.h"

int blocking_communication(algoparam_t* param, MPI_Datatype column_type, MPI_Request* reqs)
{
    // nx and ny plus ghost cells
    const int nx = param->topo_nx + 2;
    const int ny = param->topo_ny + 2;

    // Get rank again
    int rank = param->rank;

    // Create some tags for communication
    int toptag = 0;
    int downtag = 1;
    int lefttag = 2;
    int righttag = 3;

    // Calculate some indices, such that we have nice names
    int index_toprow = nx + 1;
    int index_botrow = (ny - 2) * nx + 1;

    // Destination of copy is always the ghost cells
    int index_toprow_ghost = 1;
    int index_botrow_ghost = (ny - 1) * nx + 1;

    // Only used in non-blocking communication, but the others return 0 and don't care about the loop
    int waitfor_count = 0;

    int size_row = nx - 2;

    MPI_Sendrecv(param->u + index_toprow, size_row, MPI_DOUBLE, param->topo_top, toptag, param->u + index_botrow_ghost,
                 size_row, MPI_DOUBLE, param->topo_down, toptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(param->u + index_botrow, size_row, MPI_DOUBLE, param->topo_down, downtag,
                 param->u + index_toprow_ghost, size_row, MPI_DOUBLE, param->topo_top, downtag, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

    // Exchange left most column to the left neighbors:
    MPI_Sendrecv(param->u + 1, 1, column_type, param->topo_left, lefttag, param->u + nx - 1, 1, column_type,
                 param->topo_right, lefttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Exchange right most column to the right neighbors:
    MPI_Sendrecv(param->u + nx - 2, 1, column_type, param->topo_right, righttag, param->u, 1, column_type,
                 param->topo_left, righttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return 0;
}
