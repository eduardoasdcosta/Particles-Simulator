CC=gcc
CFLAGS=-Wall -lm -fopenmp
DEPS=particle.h vectors.h grid.h
OBJ=simpar.o particle.o vectors.o grid.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simpar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm *.o simpar
