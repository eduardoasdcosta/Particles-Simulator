CC=mpicc
CFLAGS=-Wall -lm -fopenmp 
DEPS=
OBJ=simpar-mpi.o particle.o grid.o vectors.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simpar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm *.o simpar
