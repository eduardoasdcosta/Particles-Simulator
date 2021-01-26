CC=gcc
CFLAGS=-Wall -lm -fopenmp
DEPS=calculations.h particle.h grid.h
OBJ=simpar-omp.o calculations.o particle.o grid.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simpar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm *.o simpar
