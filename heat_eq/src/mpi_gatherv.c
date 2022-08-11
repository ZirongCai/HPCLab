#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "heat.h"
#include <math.h>
/*
void add_dummy_data(algoparam_t* param){
// this function is not used so far
    int lres_x = param->lresx;
    int lres_y = param->lresy;
    int res = param->visres;
    double *uvis = param->uvis;

        for (int i=0; i< res; i++){
            for (int j=0; j< res; j++){
                uvis[i*res+j] = 0.000000;
            }
    }
}
*/
void add_dummy_data(algoparam_t* param){
// this function is not used so far
    int lres_x = param->lresx;
    int lres_y = param->lresy;
    int res = param->visres;
    double *uvis = param->uvis;

    if (lres_x != res || lres_y != res){
        for (int i=lres_x; i< res; i++){
            for (int j=lres_y; j< res; j++){
                uvis[i*res+j] = 0.000000;
            }
        }    
    }
}


void print2dArray(double* array, int rows, int columns);

void all_gatherv(algoparam_t* param, MPI_Comm* comm, int size, int rank, double* gblock_ptr)
{
    // info about mpi topo
    const int nDims = 2;           
    const int dimx = param->dimx, dimy=param->dimy; 
    const int lrows  = param->visres, lcolumns  = param->visres; // lrows and lcolumns of each local block which does not contain boundaries

    // lblock is 2d data array in each local block
    double* lblock_ptr = param->uvis;

    // display local array
/*
    for (int i = 0; i < size; i++){
        if (i == rank) {
            printf("\n[Rank] of [total]: No%d of %d\n", rank, size);
            print2dArray(lblock_ptr, lrows+2, lcolumns+2);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
*/
    const int grows = lrows * dimx;          // Global array lrows
    const int gcolumns = lcolumns * dimy; // Global array lcolumns
    int gsizes[2];                     // No of elements in each dimension of the whole array
    int lsizes[2];                  // No of elements in each dimension of the subarray
    int startCoords[2];               // Starting coordinates of each subarray
    MPI_Datatype recvBlock, recvMagicBlock;

    if (rank == 0){         // For the master's eyes only
        //gblock_ptr = (double *) malloc(grows * gcolumns * sizeof(double));
        // Create a subarray (a rectangular block) datatype from a regular, 2d array
        gsizes[0] = grows;
        gsizes[1] = gcolumns;
        lsizes[0] = lrows;
        lsizes[1] = lcolumns;
        startCoords[0] = 0;
        startCoords[1] = 0;
        // At root, setup received format so that data is received in a block and placed at a right place
        MPI_Type_create_subarray(nDims, gsizes, lsizes, startCoords, MPI_ORDER_C, MPI_DOUBLE, &recvBlock);
        MPI_Type_create_resized(recvBlock, 0, lcolumns* sizeof(double), &recvMagicBlock);
        MPI_Type_commit(&recvMagicBlock);
    }

    // Format sending data so that they do not send the boundaries.
    MPI_Datatype inner_ldata_block;
    MPI_Type_vector(lrows, lcolumns, lcolumns+2, MPI_DOUBLE, &inner_ldata_block);
    MPI_Type_commit(&inner_ldata_block);

    /* MPI_Gathering.. */
    int recvCounts[size], displacements[size];

    // recvCounts: how many chunks of data each process has -- in units of blocks here --
    for (int i = 0; i < size; i++)
        recvCounts[i] = 1;

    // displacements: displacement relative to global block at which to place the
    //                             incoming data block from process i -- in block extents! --
    int index = 0;
    for (int p_row = 0; p_row < dimx; p_row++)
        for (int p_column = 0; p_column < dimy; p_column++)
            displacements[index++] = p_column  +  p_row * (lrows * dimy);

    // Gather the local arrays to the global array in the master process
    MPI_Gatherv(&lblock_ptr[lcolumns+3], 1 , inner_ldata_block, gblock_ptr, recvCounts, displacements, recvMagicBlock, 0, *comm);
    MPI_Barrier(*comm);
/*
        if (rank == 0){
        printf("Received data block at rank %d \n", rank);
        print2dArray(gblock_ptr, grows, gcolumns);
        }
  */      
        if (rank == 0) {
            MPI_Type_free(&recvMagicBlock);
            MPI_Type_free(&recvBlock);
        }
        MPI_Type_free(&inner_ldata_block);


    //return gblock_ptr;
}

void print2dArray(double* array, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++){
        for (j = 0; j < columns; j++){
            printf("%f ", array[i*columns+j]);
        }
        printf("\n");
    }
    fflush(stdout);
}
