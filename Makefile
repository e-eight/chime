CXX = g++
CXXFLAGS = -O3 -Wno-c++0x-compat -std=c++14 -fopenmp
LDLIBS = -lgomp -lm -lgsl -lgslcblas

exe = chime.exe
src = $(wildcard src/*.cpp) \
	$(wildcard src/basis/*.cpp) \
	$(wildcard src/basis/am/*.cpp) \
	$(wildcard src/basis/mcutils/*.cpp)
obj = $(src:.cpp=.o)
lib = $(wildcard *.a)

.PHONY = all tidy cleanobj cleandata

all: $(exe)

tidy:
	$(RM) $(obj) $(exe) *.out *.dat *.txt

cleanobj:
	$(RM) $(obj) $(exe)

cleandata:
	$(RM) *.dat *.txt

$(exe): $(obj)
	$(CXX) $(LDLIBS) $(CXXFLAGS) -o $@ $^
