#include "grid.h"
#include <time.h>

void print_outputs(vec2D output1, vec2D output2)
{
	printf("\n%.2lf %.2lf", output1.x, output1.y);
	printf("\n%.2lf %.2lf\n", output2.x, output2.y);
}

int main(int argc, char** argv)
{
	clock_t start = 0, end = 0;

	long seed = 0, grid_size = 0;
	long long n_par = 0, n_timesteps = 0;
	particle_t* par_array = NULL;
	long long i = 0, j = 0;
	grid** cell_grid = NULL;
	vec2D output, force;
	output.x = 0.0;
	output.y = 0.0;
	force.x = 0.0;
	force.y = 0.0;


	start = clock();
	//check if number of arguments is valid
	if(argc != 5)
	{
		fprintf(stderr, "Invalid number of arguments!\n");
		exit(1);
	}
	//get arguments from command line
	seed = atol(argv[1]);
	grid_size = atol(argv[2]);
	n_par = atoll(argv[3]);
	n_timesteps = atoll(argv[4]);

	if(grid_size < 3 || n_par < 0 || n_timesteps < 1)
	{
		fprintf(stderr, "Invalid arguments!\n");
		exit(1);
	}

	par_array = (particle_t*)malloc(n_par*sizeof(particle_t));
	if(par_array == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar.c line 44!\n");
		exit(1);
	}	

	init_particles(seed, grid_size, n_par, par_array);

	cell_grid = create_grid(grid_size);

	for(i = 0; i <= n_timesteps; i++)
	{
		#pragma omp parallel if(n_par > 1000)
		{
			update_grid(par_array, n_par, grid_size, cell_grid);

			calculate_centermass(grid_size, cell_grid);
			if(i != n_timesteps)
			{	
			    #pragma omp for private(j, force) schedule(static)
				for(j = 0; j < n_par; j++)
				{
					force = calculate_resulting_force(par_array[j], cell_grid, grid_size);

					update_particle(&par_array[j], force, grid_size);
				}
			}
		}
	}

	output = calculate_final_centermass(grid_size, cell_grid);

	printf("\n%.2lf %.2lf", par_array[0].x, par_array[0].y);

	printf("\n%.2lf %.2lf\n", output.x, output.y);

	delete_grid(grid_size, cell_grid);

	free(par_array);

	end = clock();

	fprintf(stdout, "\nExecution Time: %lf ms\n", 1000*(double) (end - start) / CLOCKS_PER_SEC);

	fflush(stdout);

	return 0;
}

