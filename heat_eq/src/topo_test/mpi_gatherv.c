#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

void print2dArray(int* array, int rows, int columns);


int main(int argc, char** argv)
{
    int size, rank;

    // Initialize MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // info about mpi topo
    int nDims = 2;           
    int dimx = 2, dimy = 2; 
    // local vars = lvars, global vars = gvars
    int lrows  = 6, lcolumns  = 8; // lrows and lcolumns of each local block


    // lblock is 2d data array in each local block
    int* lblock_ptr = (int *) malloc(lrows*lcolumns*sizeof(int)); 

    // populate arrays
    for (int i = 0; i < lrows; i++){
        for (int j = 0; j < lcolumns; j++){
            lblock_ptr[i*lcolumns+j] = i+rank;
        }
    }
    // display local array
    for (int i = 0; i < 4; i++){
        if (i == rank) {
            printf("\n[Rank] of [total]: No%d of %d\n", rank, size);
            print2dArray(lblock_ptr, lrows, lcolumns);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    int* gblock_ptr;
    //int* recvPtr;                    // Pointer to start of Global array
    int grows = (lrows-2) * dimx;          // Global array lrows
    int gcolumns = (lcolumns-2) * dimy; // Global array lcolumns
    int gsizes[2];                     // No of elements in each dimension of the whole array
    int lsizes[2];                  // No of elements in each dimension of the subarray
    int startCoords[2];               // Starting coordinates of each subarray
    MPI_Datatype recvBlock, recvMagicBlock;

    if (rank == 0){         // For the master's eyes only

        gblock_ptr = (int *) malloc(grows * gcolumns * sizeof(int));
        // Create a subarray (a rectangular block) datatype from a regular, 2d array
        gsizes[0] = grows;
        gsizes[1] = gcolumns;
        lsizes[0] = lrows-2;
        lsizes[1] = lcolumns-2;
        startCoords[0] = 0;
        startCoords[1] = 0;

        MPI_Type_create_subarray(nDims, gsizes, lsizes, startCoords, MPI_ORDER_C, MPI_INT, &recvBlock);

        // Now modify the newly created datatype to fit our needs, by specifying
        // (lower bound remains the same = 0)
        // - new extent
        // 
        MPI_Type_create_resized(recvBlock, 0, (lcolumns-2) * sizeof(int), &recvMagicBlock);

        MPI_Type_commit(&recvMagicBlock);
        //recvPtr = &gblock_ptr[0];
    }

        MPI_Datatype inner_ldata;
        MPI_Type_vector(lrows-2, lcolumns-2, lcolumns, MPI_INT, &inner_ldata);
        MPI_Type_commit(&inner_ldata);



    /* MPI_Gathering.. */
    int recvCounts[size], displacements[size];

    // recvCounts: how many chunks of data each process has -- in units of blocks here --
    for (int i = 0; i < size; i++)
        recvCounts[i] = 1;

    // dimx * dimy = np
    // displacements: displacement relative to global buffer (grow_ptr) at which to place the
    //                             incoming data block from process i -- in block extents! --
    int lcolums_last=5;
    int lrows_last=3;
    int index = 0;
    for (int p_row = 0; p_row < dimx; p_row++)
        for (int p_column = 0; p_column < dimy; p_column++)
            displacements[index++] = p_column  +  p_row * ((lrows-2) * dimy);

    // MPI_Gatherv(...) is a collective routine
    // Gather the local arrays to the global array in the master process
    // send type: MPI_INT       (a char)
    // recv type: recvMagicBlock (a block)
        MPI_Gatherv(lblock_ptr, 1, inner_ldata, //: parameters relevant to sender
                gblock_ptr, recvCounts, displacements, recvMagicBlock, 0, //: parameters relevant to receiver
                MPI_COMM_WORLD);
        //int returnvalue=mpi_gatherv(&lblock_ptr[lcolumns+1], (lrows-2) * (lcolumns-2), mpi_int, //: parameters relevant to sender
          //      gblock_ptr, recvcounts, displacements, recvmagicblock, 0, //: parameters relevant to receiver
            //    mpi_comm_world);
    // display global array
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        printf("Received data block at rank %d \n", rank);
        print2dArray(gblock_ptr, grows, gcolumns);
    }

    free(lblock_ptr);
    if (rank == 0) {
        free(gblock_ptr);
        MPI_Type_free(&recvMagicBlock);
        MPI_Type_free(&recvBlock);
    }


    MPI_Finalize();
    return 0;
}


void print2dArray(int* array, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < columns; j++){
            printf("%d ", array[i*columns+j]);
        }
        printf("\n");
    }
    fflush(stdout);
}
