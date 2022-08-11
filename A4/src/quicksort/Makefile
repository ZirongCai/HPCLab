CXX=icpc
CXXFLAGS=-Ofast -std=c++11 -march=core-avx2 -qopenmp -fno-alias -fargument-noalias-global
LDFLAGS=-qopenmp
SOURCES=quicksort.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=quicksort

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)

clean:
	rm *.o $(EXECUTABLE)
