//.h file for particle related structures and functions
#ifndef PARTICLE_H
#define PARTICLE_H

#include "vectors.h"


#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005

/* Structure containing particle information as well as their position in the grid */
typedef struct 
{
	double x,y;
    double vx,vy;
    double mass;
	int cell_line, cell_column, active;
    long long array_pos;
} particle_t;

typedef struct _proc_params
{
	int size_x, size_y, start_x, start_y, end_x, end_y, rem_x, rem_y;
    int grid_id, active;
    long long array_size;
}ProcParams;

typedef struct _particle_array
{
    particle_t *array;
    long long index;
    long long size;
}ParArray;


/* Function that creates particles (provided by the professor) */
void init_particles(long seed, long ncside, long long n_part, ParArray *par, ProcParams params);

#endif
