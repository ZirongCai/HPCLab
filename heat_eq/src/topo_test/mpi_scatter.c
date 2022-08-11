#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
//https://stackoverflow.com/questions/10788180/sending-columns-of-a-matrix-using-mpi-scatter
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
    MPI_Comm rowComm;
    MPI_Cart_sub(comm, remain_dims, &rowComm);

    //// Determind each row_rank (not 2d rank)
    int mycol;
    MPI_Comm_rank(rowComm, &mycol);

    // Create comm group by column
    int remain_dims2[2]; remain_dims2[0]=1; remain_dims2[1]=0;
    MPI_Comm colComm;
    MPI_Cart_sub(comm, remain_dims2, &colComm);

    // Determind each row_rank (not 2d rank)
    int myrow;
    MPI_Comm_rank(colComm, &myrow);

    // Block[2] == Dim;
    // square sub-block for each process 
    const int blocksize = 4;
    // size of the array need to distribute
    int globalsizes[2]={8,8};

    int *globalptr, *A; 
    int *rowptr, *rowdata;
    if (mycol == 0) {                                                                   
        if (myrow == 0){
            A=malloc(sizeof(int)*globalsizes[0]*globalsizes[1]);
            //            int A[globalsizes[0]][globalsizes[1]];
            int i,j;                                                               
            /* Initialize the matrix.  Note that C has row-major storage */        
            for (j=0; j<globalsizes[0]; j++)                                                   
                for (i=0; i<globalsizes[1]; i++)                                               
                    A[j*globalsizes[1]+i] = j;                                                   
            for (j=0; j<globalsizes[0]; j++) {                                                 
                for (i=0; i<globalsizes[1]; i++)                                               
                    printf( "%d ", A[j*globalsizes[1]+i] );                                      
                printf( "\n" );                                                    
            }                                                                      
            globalptr = &A[0];                                                    
        }                                                                          


        int sendcounts[ dim[0] ];                                                    
        int senddispls[ dim[0] ];                                                    
        senddispls[0] = 0;                                                              

        for (int row=0; row<dim[0]; row++) {                                         
            /* each processor gets blocksize rows, each of size globalsizess[1]... */    
            sendcounts[row] = blocksize*globalsizes[1];                                 
            if (row>0) 
                senddispls[row] = senddispls[row-1] + sendcounts[row-1];
        }                                                                               

        /* the last processor gets one more */                                          
        //sendcounts[dim[0]-1] += globalsizes[1];                                      


        /* allocate my rowdata, recv buf */  
        rowdata = malloc(sizeof(int)*blocksize*globalsizes[1]);
        rowptr = &rowdata[0];
        /*
           if (myrow == 0) {
           printf( "Before Output for process %d\n", rank );
           int i,j;
           for (j=0; j<blocksize; j++) {
           for (i=0; i<globalsizes[1]; i++) 
        //                    printf( "%11.11u (%d,%d) ",globalptr[globalsizes[1]*j+i],j,i );
        printf( " %d",globalptr[globalsizes[1]*j+i]);
        //printf( "%d ", recvarray[j] );
        printf( "\n" );
        }
        fflush( stdout );
        }
         */

        /* perform the scatter of rows */                                               
        MPI_Scatterv(A, sendcounts, senddispls, MPI_INT,                       
                rowdata, sendcounts[myrow], MPI_INT, 0, colComm);             

        if (rank == 2) {
            printf( "After scatter by column \n");
            printf( "Output for process %d\n", rank );
            int i,j;
            for (j=0; j<blocksize; j++) {
                for (i=0; i<globalsizes[1]; i++) 
                    //                    printf( "%11.11u (%d,%d) ",rowdata[j][i],j,i );
                    printf( " %d",rowdata[j*globalsizes[1]+i]);
                //printf( "%d ", recvarray[j] );
                printf( "\n" );
            }
            fflush( stdout );
        }
    }                                                                                   
	MPI_Barrier(colComm);
    /* Now, within each row of processors, we can scatter the columns.                
     * We can do this as we did in the previous example; create a vector              
     * (and localvector) type and scatter accordingly */                              
    //https://stackoverflow.com/questions/10788180/sending-columns-of-a-matrix-using-mpi-scatter            
    /*
       if ( isLastRow(myrow, blocks) )                                                   
       locnrows++;                                                                   
     */

    //int locnrows = blocksize;
    int locnrows = blocksize;
    int localsizes[2]={blocksize,blocksize};
    int localdata[localsizes[0]][localsizes[1]];

    MPI_Datatype vec, localvec;
	MPI_Type_vector(locnrows, 1, globalsizes[1], MPI_CHAR, &vec); 
    //MPI_Type_vector(locnrows, localsizes[1] , localsizes[1]*dim[1], MPI_INT, &vec);                     
    MPI_Type_create_resized(vec, 0, sizeof(int), &vec);                              
    MPI_Type_commit(&vec);                                                            

    MPI_Type_vector(locnrows, 1, localsizes[1], MPI_INT, &localvec);                 
    MPI_Type_create_resized(localvec, 0, sizeof(int), &localvec);                    
    MPI_Type_commit(&localvec);                                                       

    int sendcounts[ dim[1] ];                                                 
    int senddispls[ dim[1] ];                                                 
    /*
    if (mycol == 0) { 
        sendcounts[0]=localsizes[1]*blocksize;
        senddispls[0]=0;
        for (int col=1; col<dim[1]; col++) {                                  
            sendcounts[col] = localsizes[1]*blocksize;  
            senddispls[col] = col*localsizes[1];                                     
            //senddispls[col] = senddispls[col-1] + sendcounts[col-1];
        }                                                                        
    } */
 if (mycol == 0) {
        for (int col=0; col<dim[1]; col++) {
            sendcounts[col] =  blocksize;
            senddispls[col] = col*blocksize;
        }
    }
 
    rowptr = (mycol == 0) ? &(rowdata[0]) : NULL;                       

    MPI_Scatterv(rowptr, sendcounts, senddispls, MPI_INT,                            
            &(localdata[0][0]), sendcounts[mycol], MPI_INT, 0, rowComm);        


    if (rank == 0) {
        printf( "After scatter by row \n");
        printf( "Output for process %d\n", rank );
        int i,j;
        for (j=0; j<localsizes[0]; j++) {
            for (i=0; i<localsizes[1]; i++) 
                //                    printf( "%11.11u (%d,%d) ",rowdata[j][i],j,i );
                printf( " %d",localdata[j][i]);
            //printf( "%d ", recvarray[j] );
            printf( "\n" );
        }
        fflush( stdout );
    }

    // Everyone can now print their local matrix
    //MPI_Type_free(&localvec);                                                    
    //MPI_Type_free(&vec);                                                         

    MPI_Comm_free(&rowComm);                                                     
    MPI_Comm_free(&colComm);                                                     
    if (mycol==0)
        free(rowdata);

    if(rank==0){
        free(A);
    }
    MPI_Finalize();
    return 0;
}
