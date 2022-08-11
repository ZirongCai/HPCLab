K=128

# For OPT=1, GEBP_RE=1, M=32, N=32, m_c=32, k_c=32, m_r=32, n_r=16 -> 27 GFLOPS
# For OPT=1, GEBP_RE=1, M=128, N=64, m_c=32, k_c=64, m_r=32, n_r=2 -> 30 GFLOPS
#
M=32
N=32

OPT=1
# 1 for True, 0 for False
PRINT_FORMAT=1

# for OPT1
m_c=32
k_c=32
m_r=32
n_r=2

# for OPT2
BLOCK_SIZE=32

ALIGNMENT=64

CXX=icpc
CXXFLAGS+=-Ofast -std=c++11 -march=core-avx2 -fargument-noalias-global -fno-alias #-fno-alias #-fomit-frame-pointer
#CXXFLAGS+=-qopt-report-phase=vec,loop -qopt-report=5
LDFLAGS=

ifeq ($(OPT),1)
DEFINES=-DALIGNMENT=$(ALIGNMENT) -DM=$(M) -DN=$(N) -DK=$(K) -Dk_c=$(k_c) -Dm_c=$(m_c) -Dm_r=$(m_r) -Dn_r=$(n_r) -DOPT1
SOURCES=dgemm.cpp kernel1.cpp
else
DEFINES=-DALIGNMENT=$(ALIGNMENT) -DM=$(M) -DN=$(N) -DK=$(K) -DBLOCK_SIZE=$(BLOCK_SIZE)
SOURCES=dgemm.cpp kernel2.cpp
endif

ifeq ($(PRINT_FORMAT),1)
DEFINES+=-DPRINT_FORMAT
endif


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=dgemm


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@ $(DEFINES)

all: $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

clean:
	rm *.o $(EXECUTABLE)
