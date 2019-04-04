CXX = g++
CXXFLAGS = -O3 -Wno-c++0x-compat -std=c++11 -fopenmp
LDLIBS = -lgomp -lm -lgsl -lgslcblas

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
	$(CXX) -o $@ $^ $(LDLIBS) $(CXXFLAGS) 
