CXX=icpc
CXXFLAGS= -O3 -qopenmp
LDFLAGS=-qopenmp
SOURCES=pi_calculation.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=calculate_pi

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)

clean:
	rm *.o $(EXECUTABLE)
