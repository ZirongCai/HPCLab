CXX=mpicxx
CXXFLAGS+= -xCORE-AVX2 -Wall -O3 -mcmodel=medium -qopenmp -shared-intel -fno-alias -fargument-noalias-global

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
