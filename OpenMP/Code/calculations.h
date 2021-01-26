//.h file for physics related structures and functions
#ifndef CALCULATIONS_H
#define CALCULATIONS_H

#include "particle.h"
#include <math.h>

// Structure used to store positions/vectors
typedef struct 
{
	double x, y;
} vec2D;

// Structure containing the total mass,	partial x and partial y and the position of the center of mass
typedef struct 
{
	double mass_total;
	double partial_x;
	double partial_y;
	vec2D centermass_pos;

} grid;

// Function that computes the center of mass of each cell
void calculate_centermass(long grid_size, grid** cell_grid);

// Function that computes the center of mass of the overall grid
vec2D calculate_final_centermass(long grid_size, grid** cell_grid);

// Function that computes the distance between 2 points
double calculate_length(double A_x, double A_y, double B_x, double B_y);

// Function that computes the direction in the form of a vector
vec2D calculate_direction(double A_x, double A_y, double B_x, double B_y, double length);

// Function that computes the magnitude of the force
double calculate_magnitude(double massA, double massB, double distance);

// Function that computes the resulting force for each particle
vec2D calculate_resulting_force(particle_t par_array, grid** cell_grid, long grid_size);

// Function that updates the particle's information
void update_particle(particle_t *par_array, vec2D force, long grid_size);

// Function that retrieves the centers of mass of the adjacent cells
void adjust_surroundings(grid** cell_grid, long grid_size, long line, long column, double* occupancy, vec2D* centers);

#endif