#ifndef GRID_H
#define GRID_H

#include "particle.h"

#define index(i, j, grid_size) (j + (i*grid_size))

typedef struct _grid_node
{
    double centermass_x, centermass_y;
    double total_mass;
}Gnode;

void clearGrid(Gnode* grid, long grid_size);
void calculate_centermass(Gnode* grid, long grid_size, particle_t* par_array, long long n_par);
Vector2D calculate_force(particle_t par, Gnode* grid, long grid_size);
void compute_newstate(Gnode* grid, long grid_size, particle_t* par_array, long long n_par);
Vector2D calculate_finalcentermass(particle_t* par_array, long long n_par);
#endif
