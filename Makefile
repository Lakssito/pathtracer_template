CXX = g++
CXXFLAGS = -O2 -Xpreprocessor -fopenmp
OMP_INC = /opt/homebrew/opt/libomp/include
OMP_LIB = /opt/homebrew/opt/libomp/lib

all:
	$(CXX) $(CXXFLAGS) -I$(OMP_INC) -L$(OMP_LIB) -lomp -o pathtracer main.cpp && ./pathtracer
