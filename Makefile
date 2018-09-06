CC = gcc
FFLAGS = -O3 -Wall -std=c99
LFLAGS = -lm -lgsl -lgslcblas
OBJECTS = wavefunction.o main.o

main.exe: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o matelm.exe

%.o: %.c
	$(CC) $(LFLAGS) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.exe *.out *.dat
