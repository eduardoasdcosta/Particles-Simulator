//.c file for grid related functions
#ifndef GRID_H
#define GRID_H

#include "particle.h"
#include "calculations.h"


/* Function that creates 2D grid containing the particles */
grid** create_grid(long grid_size);

/* Function that updates the values of the structure of each grid */
void update_grid(particle_t* par_array, long long n_par, long grid_size, grid** cell_grid);

/* Auxiliary function that prints the information in each cell of the grid */
void print_grid(long grid_size, grid** cell_grid);

/* Function that erases the grid */
void delete_grid(long grid_size, grid** cell_grid);

#endif