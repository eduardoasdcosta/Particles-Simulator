//.c file for grid related functions
#include "grid.h"

grid** create_grid(long grid_size)
{
	int i = 0;

    grid** cell_grid = (grid**)malloc(grid_size*sizeof(grid*));

    if(cell_grid == NULL)
    {
        fprintf(stdout, "Error in create_grid function, line 12!\n");
        exit(1);
    }

    for(i = 0; i < grid_size; i++)
    {
        cell_grid[i] = (grid *)malloc(grid_size*sizeof(grid));

        if(cell_grid[i] == NULL)
	    {
	        fprintf(stdout, "Error in create_grid function, line 22!\n");
	        exit(1);
	    }
    }

    return cell_grid;

}


void update_grid(particle_t* par_array, long long n_par, long grid_size, grid** cell_grid)
{
	int i = 0, j = 0;

	// Go through every row and column 
	#pragma omp for collapse(2) private(i, j)
	for(i = 0; i < grid_size; i++)
	{
		for(j = 0; j < grid_size; j++)
		{
			// Reset values of the grid structure for each cell
    		cell_grid[i][j].mass_total = 0.0;
    		cell_grid[i][j].partial_x = 0.0;
    		cell_grid[i][j].partial_y = 0.0;
    		cell_grid[i][j].centermass_pos.x = 0.0;
    		cell_grid[i][j].centermass_pos.y = 0.0;
		}
	}
	#pragma omp for private(i) schedule(static)
	for(i = 0; i < n_par; i++)
	{
		#pragma omp atomic
			cell_grid[par_array[i].cell_line][par_array[i].cell_column].mass_total += par_array[i].m;
		#pragma omp atomic
			cell_grid[par_array[i].cell_line][par_array[i].cell_column].partial_x += par_array[i].m * par_array[i].x;
		#pragma omp atomic
			cell_grid[par_array[i].cell_line][par_array[i].cell_column].partial_y += par_array[i].m * par_array[i].y;
	}
}

void delete_grid(long grid_size, grid** cell_grid)
{
	int i = 0;

	// Free mmemory allocated to create the 2D grid
	for(i = 0; i < grid_size; i++)
	{	
		free(cell_grid[i]);
	}
	free(cell_grid);
}

