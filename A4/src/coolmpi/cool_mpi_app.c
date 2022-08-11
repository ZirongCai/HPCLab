#include <mpi.h>
//#include "/dss/dsshome1/lrz/sys/spack/release/21.1.1/opt/x86_64/intel-parallel-studio/cluster.2018.4-gcc-w3nucin/itac/2018.4.025/include/VT.h"
#include "VT.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <assert.h>

static int _event_handle;

int do_something(int work) {
  VT_begin(_event_handle); 
  usleep(work);
  VT_end(_event_handle);
}

int main(int argc, char **argv) {
  int rank, size;
  int ierr;

  int provided;
  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  assert(ierr == MPI_SUCCESS);

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(ierr == MPI_SUCCESS);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &size);
  assert(ierr == MPI_SUCCESS);
  
  char *name = "do_something"; 
  ierr = VT_funcdef(name, VT_NOCLASS, &_event_handle);
  assert(ierr == VT_OK);

  //parse arguments
  if(argc!=size+3) {
    fprintf(stderr, "Error you did not provide the correct number of arguments! Terminating...\n");
    return -1;
  }

  int dummy_buffer = -1;

  int n_its;
  n_its = atoi(argv[1]);
  
  int work_per_task;
  work_per_task = atoi(argv[2])*1000; 

  int n_tasks;
  n_tasks = atoi(argv[3+rank]);

  for(int i=0; i<n_its; i++) {
    for(int t=0; t<n_tasks; t++) {
      if(rank==1)
        do_something(2*work_per_task);
      else
        do_something(work_per_task);
    }
    if(rank==0) dummy_buffer = i;
    MPI_Allreduce( &i, &dummy_buffer, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
  }

  MPI_Finalize();
}
