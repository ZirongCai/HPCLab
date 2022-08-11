CXX = mpicxx
CXX_FLAGS = --std=c++17 -Wall -Wextra -mavx -O3 -g -DOMPI_SKIP_MPICXX
#LDFLAGS=-qopenmp
SOURCES=broadcast.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=broadcast

# this compiler definition is needed to silence warnings caused by the openmpi CXX
# bindings that are deprecated. This is needed on gcc 8 forward.
# see: https://github.com/open-mpi/ompi/issues/5157

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX)  -c $(CXXFLAGS) $< -o $@ $(DEFINES)

clean:
	rm *.o $(EXECUTABLE)

