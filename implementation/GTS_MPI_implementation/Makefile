# Include files
SOURCES=main.c computeFlux.c computeRoeSolver.c updateFlux.c visual.c helper.c

# Compiler
# --------
CC=mpicc

CFLAGS=-Wall -std=c99#-Werror -pedantic

# Linker flags
# ------------
LDFLAGS=-lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=swe_sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) *.vtk


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

