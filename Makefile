CC = gcc
FFLAGS = -O3 -Wall -std=c99 -fopenmp
LFLAGS = -lgomp -lm -lgsl -lgslcblas -lcuba
OBJECTS = tdho.o constants.o utility.o main.o \
	operator_r_sq.o 
SHAREDLIBS = libcuba.a

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) $(SHAREDLIBS) -o matelm.exe

%.o: %.c
	$(CC) $(LFLAGS) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.out *.dat
