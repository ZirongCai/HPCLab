#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>

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
    MPI_Cart_coords(comm, rank, 2, coord);

    // Create comm group by row
    int remain_dims[2]; remain_dims[0]=0; remain_dims[1]=1;
    MPI_Comm comm_row;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    //// Determind each row_rank (not 2d rank)
    int my_row_rank;
    MPI_Comm_rank(comm_row, &my_row_rank);
    int row_coord;
    MPI_Cart_coords(comm_row, my_row_rank, 1, &row_coord);

    // Create comm group by column
    int remain_dims2[2]; remain_dims2[0]=1; remain_dims2[1]=0;
    MPI_Comm comm_column;
    MPI_Cart_sub(comm, remain_dims2, &comm_column);

    // Determind each row_rank (not 2d rank)
    int my_column_rank;
    MPI_Comm_rank(comm_column, &my_column_rank);
    int column_coord;
    MPI_Cart_coords(comm_column, my_column_rank, 1, &column_coord);
/*
    if (rank == 1)
    {
        printf("Rank %d coordinates are %d %d\n",  rank, coord[0], coord[1]);
        printf("Row_rank %d coordinates are %d\n", my_row_rank, row_coord);
    }
*/
    // Test with Scatter
    if (coord[0]==0){
        int nx, ny;
        nx = 4; ny=8;
        int gsize = 2;
        int root, *sendbuf, myrank, bufsize, *stride; 
        int recvarray[2][ny]; 
        MPI_Datatype rtype; 
        int i, *displs, *scounts, offset;

        stride = (int *)malloc(gsize*sizeof(int));
        displs = (int *)malloc(gsize*sizeof(int)); 
        scounts = (int *)malloc(gsize*sizeof(int)); 
        offset = 0; 


        if (my_row_rank == 0) {
            int A[nx][ny];
            int i,j;
            /* Initialize the matrix.  Note that C has row-major storage */
            for (j=0; j<nx; j++) 
                for (i=0; i<ny; i++)
                    A[j][i] = i;
            for (j=0; j<nx; j++) {
                for (i=0; i<ny; i++) 
                    printf( "%d ", A[j][i] );
                printf( "\n" );
            } 
            sendbuf = &A[0][0];
        }

        //stride={0,8,8,8};
        for (i=0; i<gsize; ++i) { 
            displs[i] = offset; 
            //offset += stride[i]; 
            offset += 2*ny; 
            scounts[i] = 2*ny; //nx 
        } 
        /* Create datatype for the column we are receiving 
         */ 
        /* 
        //MPI_Type_vector( 2, 4, 4, MPI_INT, &rtype); 
        int mat_size[2]={nx,ny};
        int submat_size[2]={2,ny};
        int submat_start[2]={rank*2,0};
        MPI_Type_create_subarray(2, mat_size, submat_size,submat_start, MPI_ORDER_C, MPI_INT, &rtype);
        MPI_Type_commit( &rtype ); 
        //rptr = &recvarray[0][rank]; 
        MPI_Scatterv( sendbuf, scounts, displs, MPI_INT, &recvarray[0][0], 8, rtype, 0, MPI_COMM_WORLD);
         */
        MPI_Scatterv(sendbuf, scounts, displs, MPI_INT, recvarray, 16, MPI_INT, 0, comm_row);
        /* Everyone can now print their local matrix */

        if (my_row_rank == 1) {
            printf( "Output for process %d\n", rank );
            int i,j;
            for (j=0; j<2; j++) {
                for (i=0; i<ny; i++) 
                    printf( "%d ", recvarray[j][i] );
                //printf( "%d ", recvarray[j] );
                printf( "\n" );
            }
            fflush( stdout );
        }
    }
    MPI_Finalize();
    return 0;
}
