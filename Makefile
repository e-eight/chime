CXX = g++
FFLAGS = -O3 -Wall -std=c++14 -fopenmp
LFLAGS = -lgomp -lm -lgsl -lgslcblas

exe = chime.exe
src = $(wildcard src/*.cpp) \
	$(wildcard src/basis/*.cpp) \
	$(wildcard src/basis/am/*.cpp) \
	$(wildcard src/basis/mcutils/*.cpp)
obj = $(src:.cpp=.o)
lib = $(wildcard *.a)

.PHONY = all clean

all: $(exe)

clean:
	rm -rf $(obj) $(exe) *.out *.dat

$(exe): $(obj)
	$(CXX) -o $@ $^ $(LFLAGS) $(FFLAGS)
