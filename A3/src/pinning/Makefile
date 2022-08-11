CXX=mpicxx
CXXFLAGS+=-qopenmp
LDFLAGS=-qopenmp

SOURCES=main.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=pinning


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)

all: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

clean:
	rm *.o $(EXECUTABLE)
