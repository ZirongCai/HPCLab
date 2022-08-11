# Adjust these values, if necessary
M=32
N=2
K=64
MC=32
ALIGNMENT=64

TASK=1 # comment it out if you want to have Loop-based

CXX=icpc
CXXFLAGS=-Ofast -std=c++11 -march=core-avx2 -qopenmp -fno-alias -fargument-noalias-global
#CXXFLAGS+=-qopt-report-phase=vec,loop -qopt-report=5

CXXFLAGS+=-DNOFT # No First-Touch, Comment it out if you want to have First-Touch

LDFLAGS=-qopenmp #-lxsmm -lblas
DEFINES=-DALIGNMENT=$(ALIGNMENT) -DM=$(M) -DN=$(N) -DK=$(K) -DMC=$(MC)
ifdef TASK
SOURCES=dgemm.cpp kernel_task.cpp
else
SOURCES=dgemm.cpp kernel.cpp
endif

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=dgemm


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)
	
clean:
	rm *.o $(EXECUTABLE)

