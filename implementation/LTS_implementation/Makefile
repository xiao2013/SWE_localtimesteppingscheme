# Include files
SOURCES=main.c computeFlux.c computeTimeStep.c computeRoeSolver.c updateFlux.c visual.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-Wall #-Werror -pedantic

# Linker flags
# ------------
LDFLAGS= 

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=swe_sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

