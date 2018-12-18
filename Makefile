CXX = g++
FFLAGS = -O3 -Wall -std=c++14 -fopenmp
LFLAGS = -lgomp -lm -lgsl -lgslcblas

exe = ndrel
src = $(wildcard **/*.cpp) 
obj = $(src:.cpp=.o)
lib = $(wildcard *.a)

.PHONY = all clean

all: $(exe)

clean:
	rm -rf $(obj) *.exe *.out *.dat

$(exe): $(obj)
	$(CXX) -o $@ S^ $(LFLAGS) $(FFLAGS) 
