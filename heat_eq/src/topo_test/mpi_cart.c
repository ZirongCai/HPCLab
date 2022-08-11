#include<mpi.h>
#include<stdio.h>

/* A two-dimensional processes in a 2x2 grid */
int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Comm comm;
    int dim[2], period[2], reorder;
    int coord[2], id;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    dim[0]=2; dim[1]=2;
    period[0]=0; period[1]=0;
    reorder=0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);
    int cart_rank;
    MPI_Comm_rank(comm, & cart_rank);

    // Create comm group by row
    int remain_dims[2]; remain_dims[0]=0; remain_dims[1]=1;
    MPI_Comm comm_row;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    // Determind each row_rank (not 2d rank)
    int my_row_rank;
    MPI_Comm_rank(comm_row, &my_row_rank);
    MPI_Cart_coords(comm, rank, 2, coord);

    if (rank == 2)
    {
        printf("Rank %d coordinates are %d %d\n",  rank, coord[0], coord[1]);
        printf("cart_rank %d \n",  cart_rank);
        printf("Row_rank %d \n", my_row_rank);
        //MPI_Cart_coords(comm_row, rank, 2, coord);
        //printf("Rank %d coordinates are %d %d\n", rank, coord[0], coord[1]);
    }

    MPI_Finalize();
    return 0;
}
