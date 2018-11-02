CC = gcc
FFLAGS = -O3 -Wall -std=c99 -fopenmp
LFLAGS = -lgomp -lm -lgsl -lgslcblas -lcuba

MKDIR_P = mkdir -p

SOURCEDIR = src
BUILDDIR = build
LIBDIR = local/lib

EXECUTABLE = matelm.exe
SOURCES = $(wildcard $(SOURCEDIR)/*.c)
OBJECTS = $(patsubst $(SOURCEDIR)/%.c, $(BUILDDIR)/%.o, $(SOURCES)) 
SHAREDLIBS = $(wildcard $(LIBDIR)/*.a)

.PHONY = directories all clean

all: $(EXECUTABLE)

directories: $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR) *.exe *.out *.dat

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(SHAREDLIBS) $(LFLAGS) -o $@

$(OBJECTS): $(BUILDDIR)/%.o: $(SOURCEDIR)/%.c
	$(CC) $(LFLAGS) $(FFLAGS) -c $< -o $@

$(BUILDDIR):
	$(MKDIR_P) $(BUILDDIR)
