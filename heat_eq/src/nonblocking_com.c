
#include "communication.h"

int nonblocking_communication(algoparam_t* param, MPI_Datatype column_type, MPI_Request* reqs)
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

    // IMPORTANT!! Read first!
    // All sends, send non-ghost cells!
    // All recvs, save into the ghost cells
    // This is always the case, so some rank sends it's non ghost cells to the neighbours
    // and receives the (non-ghost) cells from the neighbours, which the rank stores into
    // its own ghost cells and starts communication again

    MPI_Request req = MPI_REQUEST_NULL;

    // Exchange the top row with top neighbors
    MPI_Isend(param->u + index_toprow, size_row, MPI_DOUBLE, param->topo_top, toptag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;
    MPI_Irecv(param->u + index_botrow_ghost, size_row, MPI_DOUBLE, param->topo_down, toptag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;

    // Exchange the BOTTOM row with bottom neighbors
    MPI_Isend(param->u + index_botrow, size_row, MPI_DOUBLE, param->topo_down, downtag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;
    MPI_Irecv(param->u + index_toprow_ghost, size_row, MPI_DOUBLE, param->topo_top, downtag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;

    // Exchange LEFT most column to the left neighbors:
    MPI_Isend(param->u + 1, 1, column_type, param->topo_left, lefttag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;
    MPI_Irecv(param->u + nx - 1, 1, column_type, param->topo_right, lefttag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;

    // Exchange RIGHT most column to the right neighbors:
    MPI_Isend(param->u + nx - 2, 1, column_type, param->topo_right, righttag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;
    MPI_Irecv(param->u, 1, column_type, param->topo_left, righttag, MPI_COMM_WORLD, &req);
    reqs[waitfor_count++] = req;

    return waitfor_count;
}
