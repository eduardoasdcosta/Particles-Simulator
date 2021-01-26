//.c file for physics related functions
#include "calculations.h"


void calculate_centermass(long grid_size, grid** cell_grid)
{
	int i = 0, j = 0;

	// Go through every row and column
	#pragma omp for collapse(2) private(i,j) schedule(static)
	for(i = 0; i < grid_size; i++)
	{
		for(j = 0; j < grid_size; j++)
		{
			// If there are no particles in the cell, move to the next one
			if(cell_grid[i][j].mass_total == 0)
			{
				continue;
			}
			else // Calculate the position of center of mass by dividing the partial x and y summations with the total mass of the cell
			{	
				cell_grid[i][j].centermass_pos.x = cell_grid[i][j].partial_x / cell_grid[i][j].mass_total;
				cell_grid[i][j].centermass_pos.y = cell_grid[i][j].partial_y / cell_grid[i][j].mass_total;
			}
		}
	}
}

vec2D calculate_final_centermass(long grid_size, grid** cell_grid)
{

	int i = 0,j = 0;
	double total_mass = 0.0;
	vec2D centermass_pos;

	centermass_pos.x = 0.0;
	centermass_pos.y = 0.0;

	// Go through every row and column
	for(i = 0; i < grid_size; i++)
	{
		for(j = 0; j < grid_size; j++)
		{

			// If there are no particles in the cell, move to the next one
			if(cell_grid[i][j].mass_total == 0)
			{
				continue;
			}
			else 
			{	
				// Add mass of every filled cell 
				total_mass += cell_grid[i][j].mass_total;

				// Calculate partial x and y summations 
				centermass_pos.x += cell_grid[i][j].mass_total * cell_grid[i][j].centermass_pos.x;
				centermass_pos.y += cell_grid[i][j].mass_total * cell_grid[i][j].centermass_pos.y;
			}
		}
	}

	// Calculate the position of center of mass by dividing the partial x and y summations with the total mass of the grid
	centermass_pos.x = centermass_pos.x / total_mass;
	centermass_pos.y = centermass_pos.y / total_mass;

	return(centermass_pos);

}

double calculate_length(double A_x, double A_y, double B_x, double B_y)
{
	return sqrt(pow(A_x - B_x, 2) + pow(A_y - B_y, 2));
}

vec2D calculate_direction(double A_x, double A_y, double B_x, double B_y, double length)
{
	vec2D direction;

	direction.x = (B_x - A_x) / length;
	direction.y = (B_y - A_y) / length;

	return direction;
}

double calculate_magnitude(double massA, double massB, double distance)
{

	return (G * (massA * massB) / pow(distance, 2));
}


vec2D calculate_resulting_force(particle_t par_array, grid** cell_grid, long grid_size)
{
	int i = 0, j = 0;
	double length = 0.0, magnitude = 0.0;
	long line = 0, column = 0;
	vec2D centers[9];
	double occupancy[9];
	vec2D force;
	vec2D direction;
	force.x = 0.0;
	force.y = 0.0;
	direction.x = 0.0;
	direction.y = 0.0;

	for(i = 0; i < 9; i++)
	{
		occupancy[i] = 0.0;
		centers[i].x = 0.0;
		centers[i].y = 0.0;
	}

	// Store the particle's location in the grid
	line = par_array.cell_line;
	column = par_array.cell_column;

	// Store the position and mass of the center of mass in that cell, respectively
	centers[4] = cell_grid[line][column].centermass_pos;
	occupancy[4] = cell_grid[line][column].mass_total;

	// Call function that retrieves the centers of mass of the adjacent cells, considering cells "outside" de grid
	adjust_surroundings(cell_grid, grid_size, line, column, occupancy, centers);

	// Go through every adjacent cell and including own cell
	for(j = 0; j < 9; j++)
	{	
		if(occupancy[j] == 0)
		{	
			// If the adjacent cell is not occupied by any particle, move to next adjacent cell
			continue;
		}
		else
		{
			// If it is occupied, calculate length to the center of mass
			length = calculate_length(par_array.x, par_array.y, centers[j].x, centers[j].y);

			// If less than predefined threshold, move to next adjacent cell
			if(length < EPSLON)
			{
				continue;
			}
			else
			{
				// Compute direction of motion
				direction = calculate_direction(par_array.x, par_array.y, centers[j].x, centers[j].y, length);

				// Compute force magnitude
				magnitude = calculate_magnitude(par_array.m, occupancy[j], length);

				// Compute force vector
				force.x += magnitude * direction.x;
				force.y += magnitude * direction.y;
			}
		}
	}

	return force;
}


void adjust_surroundings(grid** cell_grid, long grid_size, long line, long column, double occupancy[], vec2D centers[])
{
	// Cell on the inside layers
	if((line <= grid_size - 2 && line >= 1) && (column <= grid_size - 2 && column >= 1))
	{
		centers[0] = cell_grid[line-1][column-1].centermass_pos;
		occupancy[0] = cell_grid[line-1][column-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][column+1].centermass_pos;
		occupancy[2] = cell_grid[line-1][column+1].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;

		centers[6] = cell_grid[line+1][column-1].centermass_pos;
		occupancy[6] = cell_grid[line+1][column-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][column+1].centermass_pos;
		occupancy[8] = cell_grid[line+1][column+1].mass_total;
	}
	// Cell in the top left corner
	else if(line + 1 >= grid_size && column - 1 < 0) 							// top left corner
	{
		centers[0] = cell_grid[line-1][grid_size-1].centermass_pos;
		centers[0].x--;
		occupancy[0] = cell_grid[line-1][grid_size-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][column+1].centermass_pos;
		occupancy[2] = cell_grid[line-1][column+1].mass_total;

		centers[3] = cell_grid[line][grid_size-1].centermass_pos;
		centers[3].x--;
		occupancy[3] = cell_grid[line][grid_size-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;
	
		centers[6] = cell_grid[0][grid_size-1].centermass_pos;
		centers[6].x--;
		centers[6].y++;
		occupancy[6] = cell_grid[0][grid_size-1].mass_total;

		centers[7] = cell_grid[0][column].centermass_pos;
		centers[7].y++;
		occupancy[7] = cell_grid[0][column].mass_total;

		centers[8] = cell_grid[0][column+1].centermass_pos;
		centers[8].y++;
		occupancy[8] = cell_grid[0][column+1].mass_total;
	}		
	// Cell in the top edge 
	else if(line+1 >= grid_size && column+1 < grid_size && column-1 >=0)    
	{
		centers[0] = cell_grid[line-1][column-1].centermass_pos;
		occupancy[0] = cell_grid[line-1][column-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][column+1].centermass_pos;
		occupancy[2] = cell_grid[line-1][column+1].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;

		centers[6] = cell_grid[0][column-1].centermass_pos;
		centers[6].y++;
		occupancy[6] = cell_grid[0][column-1].mass_total;

		centers[7] = cell_grid[0][column].centermass_pos;
		centers[7].y++;
		occupancy[7] = cell_grid[0][column].mass_total;

		centers[8] = cell_grid[0][column+1].centermass_pos;
		centers[8].y++;
		occupancy[8] = cell_grid[0][column+1].mass_total;
	}
	// Cell in the top right corner
	else if(line+1 >= grid_size && column+1 >= grid_size)		
	{
		centers[0] = cell_grid[line-1][column-1].centermass_pos;
		occupancy[0] = cell_grid[line-1][column-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][0].centermass_pos;
		centers[2].x++;
		occupancy[2] = cell_grid[line-1][0].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][0].centermass_pos;
		centers[5].x++;
		occupancy[5] = cell_grid[line][0].mass_total;

		centers[6] = cell_grid[0][column-1].centermass_pos;
		centers[6].y++;
		occupancy[6] = cell_grid[0][column-1].mass_total;

		centers[7] = cell_grid[0][column].centermass_pos;
		centers[7].y++;
		occupancy[7] = cell_grid[0][column].mass_total;

		centers[8] = cell_grid[0][0].centermass_pos;
		centers[8].x++;
		centers[8].y++;
		occupancy[8] = cell_grid[0][0].mass_total;
	}
	// Cell on the left edge
	else if(column-1 < 0 && line-1 >=0 && line+1 < grid_size)			
	{
		centers[0] = cell_grid[line-1][grid_size-1].centermass_pos;
		centers[0].x--;
		occupancy[0] = cell_grid[line-1][grid_size-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][column+1].centermass_pos;
		occupancy[2] = cell_grid[line-1][column+1].mass_total;

		centers[3] = cell_grid[line][grid_size-1].centermass_pos;
		centers[3].x--;
		occupancy[3] = cell_grid[line][grid_size-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;

		centers[6] = cell_grid[line +1][grid_size-1].centermass_pos;
		centers[6].x--;
		occupancy[6] = cell_grid[line+1][grid_size-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][column+1].centermass_pos;
		occupancy[8] = cell_grid[line+1][column+1].mass_total;
	}
	// Cell on the right edge
	else if(column+1 >= grid_size && line+1 < grid_size && line-1 >=0) 		
	{
		centers[0] = cell_grid[line-1][column-1].centermass_pos;
		occupancy[0] = cell_grid[line-1][column-1].mass_total;

		centers[1] = cell_grid[line-1][column].centermass_pos;
		occupancy[1] = cell_grid[line-1][column].mass_total;

		centers[2] = cell_grid[line-1][0].centermass_pos;
		centers[2].x++;
		occupancy[2] = cell_grid[line-1][0].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][0].centermass_pos;
		centers[5].x++;
		occupancy[5] = cell_grid[line][0].mass_total;

		centers[6] = cell_grid[line+1][column-1].centermass_pos;
		occupancy[6] = cell_grid[line+1][column-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][0].centermass_pos;
		centers[8].x++;
		occupancy[8] = cell_grid[line+1][0].mass_total;
	}
	// Cell on the bottom left corner
	else if(line-1 < 0 && column-1 <0)										
	{
		centers[0] = cell_grid[grid_size-1][grid_size-1].centermass_pos;
		centers[0].x--;
		centers[0].y--;
		occupancy[0] = cell_grid[grid_size-1][grid_size-1].mass_total;

		centers[1] = cell_grid[grid_size -1][column].centermass_pos;
		centers[1].y--;
		occupancy[1] = cell_grid[grid_size-1][column].mass_total;

		centers[2] = cell_grid[grid_size-1][column+1].centermass_pos;
		centers[2].y--;
		occupancy[2] = cell_grid[grid_size-1][column+1].mass_total;

		centers[3] = cell_grid[line][grid_size-1].centermass_pos;
		centers[3].x--;
		occupancy[3] = cell_grid[line][grid_size-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;

		centers[6] = cell_grid[line+1][grid_size-1].centermass_pos;
		centers[6].x--;
		occupancy[6] = cell_grid[line+1][grid_size-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][column+1].centermass_pos;
		occupancy[8] = cell_grid[line+1][column+1].mass_total;
	}
	// Cell on the bottom edge
	else if(line-1 < 0 && column-1 >=0 && column+1 < grid_size) 			
	{
		centers[0] = cell_grid[grid_size-1][column-1].centermass_pos;
		centers[0].y--;
		occupancy[0] = cell_grid[grid_size-1][column-1].mass_total;

		centers[1] = cell_grid[grid_size -1][column].centermass_pos;
		centers[1].y--;
		occupancy[1] = cell_grid[grid_size-1][column].mass_total;

		centers[2] = cell_grid[grid_size-1][column+1].centermass_pos;
		centers[2].y--;
		occupancy[2] = cell_grid[grid_size-1][column+1].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][column+1].centermass_pos;
		occupancy[5] = cell_grid[line][column+1].mass_total;

		centers[6] = cell_grid[line+1][column-1].centermass_pos;
		occupancy[6] = cell_grid[line+1][column-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][column+1].centermass_pos;
		occupancy[8] = cell_grid[line+1][column+1].mass_total;
	}
	// Cell on the bottom right corner
	else if(line-1 < 0 && column+1 >= grid_size)							
	{
		centers[0] = cell_grid[grid_size-1][column-1].centermass_pos;
		centers[0].y--;
		occupancy[0] = cell_grid[grid_size-1][column-1].mass_total;

		centers[1] = cell_grid[grid_size-1][column].centermass_pos;
		centers[1].y--;
		occupancy[1] = cell_grid[grid_size-1][column].mass_total;

		centers[2] = cell_grid[grid_size-1][0].centermass_pos;
		centers[2].x++;
		centers[2].y--;
		occupancy[2] = cell_grid[grid_size-1][0].mass_total;

		centers[3] = cell_grid[line][column-1].centermass_pos;
		occupancy[3] = cell_grid[line][column-1].mass_total;

		centers[5] = cell_grid[line][0].centermass_pos;
		centers[5].x++;
		occupancy[5] = cell_grid[line][0].mass_total;

		centers[6] = cell_grid[line+1][column-1].centermass_pos;
		occupancy[6] = cell_grid[line+1][column-1].mass_total;

		centers[7] = cell_grid[line+1][column].centermass_pos;
		occupancy[7] = cell_grid[line+1][column].mass_total;

		centers[8] = cell_grid[line+1][0].centermass_pos;
		centers[8].x++;
		occupancy[8] = cell_grid[line+1][0].mass_total;
	}
	// In case of error
	else																	
	{
		fprintf(stdout, "Error in function adjust_surroundings!\n");
		exit(1);
	}
}

void update_particle(particle_t *particle, vec2D force, long grid_size)
{
	vec2D acceleration;

	// Compute acceleration vector from force vector and particle's mass
	acceleration.x = force.x/particle->m;
	acceleration.y = force.y/particle->m;

	// Compute and update the particle's velocity vector
	particle->vx = particle->vx + acceleration.x;
	particle->vy = particle->vy + acceleration.y;

	// Compute and update the particle's position
	particle->x = particle->x + particle->vx + 0.5*acceleration.x;
	particle->y = particle->y + particle->vy + 0.5*acceleration.y;

	// Check if particle's new position is outside the grid, if so, correct position 
	if(particle->x < 0)
		particle->x++;
	else if(particle->x > 1)
		particle->x--;
		
	if(particle->y < 0)
		particle->y++;
	else if(particle->y > 1)
		particle->y--;

	particle->cell_line = (particle->y) * grid_size;
	particle->cell_column = (particle->x) * grid_size;

}


