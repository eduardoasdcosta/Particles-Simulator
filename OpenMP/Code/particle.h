//.h file for particle related structures and functions
#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <stdio.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

/* Structure containing particle information as well as their position in the grid */
typedef struct 
{
	double x, y, vx, vy, m;
	int cell_line, cell_column;
} particle_t;

/* Function that creates particles (provided by the professor) */
void init_particles(long seed, long ncside, long long n_part, particle_t *par);

#endif