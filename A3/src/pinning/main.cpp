#include <iostream>
#include <string>
#include <cstdlib>
#include <sched.h>
#include <mpi.h>
#include <omp.h>

#define MPI_CHECK check_mpi_status(status, __FILE__, __LINE__)
void check_mpi_status(int status, std::string file_name, int line_num) {
    if (status != MPI_SUCCESS) {
        std::cout << "MPI::ERROR::(" << status << "):" 
                  << file_name << ":" << line_num
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}
int status = 0;

int main(int argc, char *argv[]) {

    status = MPI_Init (&argc, &argv); MPI_CHECK;

    int rank = -1, size = -1;

    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_CHECK;
    status = MPI_Comm_size(MPI_COMM_WORLD, &size); MPI_CHECK;

    int left = -1;
    if (rank > 0) {
        status = MPI_Recv(&left, 
                          1, 
                          MPI_INT, 
                          rank - 1, 
                          0,
                          MPI_COMM_WORLD, 
                          MPI_STATUS_IGNORE); MPI_CHECK;
    }
   
    unsigned flag = 0;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        while(thread_id != flag) {
            #pragma omp flush(flag)
        }

        std::cout << "rank=" << rank << "; " 
                  << "thread=" << thread_id << "; " 
                  << "cpu id=" << sched_getcpu() 
                  << std::endl;
        ++flag;
    }

    if (rank < (size - 1)) {
        status =  MPI_Send(&rank, 
                           1, 
                           MPI_INT, 
                           rank + 1, 
                           0, 
                           MPI_COMM_WORLD); MPI_CHECK;
    }


    status = MPI_Finalize(); MPI_CHECK;
    return 0;
}
