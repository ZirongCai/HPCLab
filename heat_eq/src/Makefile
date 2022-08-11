# CFLAGS += -O3 $(PAPI_INC) -std=c11 -fopenmp
CFLAGS += -O3 $(PAPI_INC) -std=c11 -O3 -xCORE-AVX512 -fno-alias -qopt-zmm-usage=high -funroll-all-loops -qopenmp -qopt-streaming-stores=always


MPICC = mpicc

all: heat 

heat : heat.o input.o misc.o timing.o split_topology.o jacobi_iters.o relax_jacobi.o blocking_com.o nonblocking_com.o mpi_gatherv.o
	$(MPICC) ${INC} ${MPI_INC} $(CFLAGS) -o heat $+ -lm  $(PAPI_LIB)

%.o : %.c %.h
	$(MPICC) ${INC} ${MPI_INC} $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o heat *~ *.ppm

remake : clean all
