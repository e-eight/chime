CC = gcc
FFLAGS = -O3 -Wall -std=c99 -fopenmp
LFLAGS = -lgomp -lm -lgsl -lgslcblas -lcuba
OBJECTS = libcuba.a tdho.o main.o utility.o operator_id.o

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o matelm.exe

%.o: %.c
	$(CC) $(LFLAGS) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.out *.dat
