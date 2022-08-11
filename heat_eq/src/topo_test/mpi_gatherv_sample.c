#include <stdio.h>
  #include <stdlib.h>
  #include <mpi.h>
  #include <time.h>

  void print2dCharArray(int** array, int rows, int columns);


  int main(int argc, char** argv)
  {
    int master = 0, np, rank;
    char version[10];
    char processorName[20];
    int strLen[10];

    // Initialize MPI environment
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    //if (np != 12) { MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // We need a different seed for each process
    //srand(time(0) ^ (rank * 33 / 4));

    int nDims = 2;               // array dimensions
    int rows  = 4, columns  = 6; // rows and columns of each block
    int prows = 2, pcolumns = 2; // rows and columns of blocks. Each block is handled by 1 process

    int* pre_grid = (int*) malloc(rows * columns * sizeof(int));
    int** grid = (int**) malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++)
      grid[i] = &(pre_grid[i * columns]);

    int** universe = NULL;           // Global array
    int* pre_universe = NULL;
    int* recvPtr;                    // Pointer to start of Global array
    int Rows = rows * prows;          // Global array rows
    int Columns = columns * pcolumns; // Global array columns
    int sizes[2];                     // No of elements in each dimension of the whole array
    int subSizes[2];                  // No of elements in each dimension of the subarray
    int startCoords[2];               // Starting coordinates of each subarray
    MPI_Datatype recvBlock, recvMagicBlock;

    if (rank == master){         // For the master's eyes only
  /*    universe = malloc(Rows * sizeof(int*));*/
  /*    for (int i = 0; i < Rows; i++)*/
  /*      universe[i] = malloc(Columns * sizeof(int));*/

      pre_universe = (int*) malloc(Rows * Columns * sizeof(int));
      universe = (int**) malloc(Rows * sizeof(int*));
      for (int i = 0; i < Rows; i++) {
          universe[i] = &(pre_universe[i * Columns]);
      }



      // Create a subarray (a rectangular block) datatype from a regular, 2d array
      sizes[0] = Rows;
      sizes[1] = Columns;
      subSizes[0] = rows;
      subSizes[1] = columns;
      startCoords[0] = 0;
      startCoords[1] = 0;

      MPI_Type_create_subarray(nDims, sizes, subSizes, startCoords, MPI_ORDER_C, MPI_INT, &recvBlock);

      // Now modify the newly created datatype to fit our needs, by specifying
      // (lower bound remains the same = 0)
      // - new extent
      // The new region / block will now "change" sooner, as soon as we reach a region of elements
      //         occupied by a new block, ie. every: (columns) * sizeof(elementType) =
      MPI_Type_create_resized(recvBlock, 0, columns * sizeof(int), &recvMagicBlock);

      MPI_Type_commit(&recvMagicBlock);
      recvPtr = &universe[0][0];
    }

    // populate arrays
    for (int y = 0; y < rows; y++){
      for (int x = 0; x < columns; x++){
        grid[y][x] = rank;
      }
    }


    // display local array
    for (int i = 0; i < np; i++){
      if (i == rank) {
        printf("\n[Rank] of [total]: No%d of %d\n", rank, np);
        print2dCharArray(grid, rows, columns);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }


    /* MPI_Gathering.. */
    int recvCounts[np], displacements[np];

    // recvCounts: how many chunks of data each process has -- in units of blocks here --
    for (int i = 0; i < np; i++)
      recvCounts[i] = 1;

    // prows * pcolumns = np
    // displacements: displacement relative to global buffer (universe) at which to place the
    //                             incoming data block from process i -- in block extents! --
    int index = 0;
    for (int p_row = 0; p_row < prows; p_row++)
      for (int p_column = 0; p_column < pcolumns; p_column++)
        displacements[index++] = p_column  +  p_row * (rows * pcolumns);

    // MPI_Gatherv(...) is a collective routine
    // Gather the local arrays to the global array in the master process
    // send type: MPI_INT       (a int)
    // recv type: recvMagicBlock (a block)
    MPI_Gatherv(&grid[0][0], rows * columns, MPI_CHAR, //: parameters relevant to sender
            recvPtr, recvCounts, displacements, recvMagicBlock, master, //: parameters relevant to receiver
            MPI_COMM_WORLD);

    // display global array
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == master){
      printf("\n---Global Array---\n");
      print2dCharArray(universe, Rows, Columns);
    }

    free(grid[0]);
    free(grid);
    if (rank == master) {
      free(universe[0]);
      free(universe);
      MPI_Type_free(&recvMagicBlock);
      MPI_Type_free(&recvBlock);
    }


    MPI_Finalize();
    return 0;
  }


  void print2dCharArray(int** array, int rows, int columns)
  {
    int i, j;
    for (i = 0; i < rows; i++){
      for (j = 0; j < columns; j++){
        printf("%d ", array[i][j]);
      }
      printf("\n");
    }
    fflush(stdout);
  }
