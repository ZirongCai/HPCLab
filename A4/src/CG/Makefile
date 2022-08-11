#CXX=mpicxx
#CXX=scalasca -instrument --mpp=mpi mpicxx # For scalasca automatically
#CXX=scalasca -instrument --nocompiler --mpp=mpi mpicxx # For scalasca manually
CXX=scalasca -instrument --nocompiler --thread=omp --mpp=mpi mpicxx # For hybrid

CXXFLAGS+= -xCORE-AVX2 -Wall -Ofast -mcmodel=medium -qopenmp -shared-intel -fno-alias -fargument-noalias-global 
#CXXFLAGS+=-qopt-report-phase=vec,loop -qopt-report=5

LDFLAGS= -qopenmp

SOURCES=parallel_solver.cpp general.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=parallel

SOURCES2=serial_code.cpp
OBJECTS2=$(SOURCES2:.cpp=.o)
EXECUTABLE2=serial

%.o: %.cpp 
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)

all: $(OBJECTS) $(OBJECTS2)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)
	$(CXX) $(LDFLAGS) $(OBJECTS2) -o $(EXECUTABLE2)

cleandata:
	rm *.gnuplot
clean:
	rm *.o $(EXECUTABLE) $(EXECUTABLE2)
