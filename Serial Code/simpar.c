#include "grid.h"

int main(int argc, char** argv)
{
    int n = 0;
	long seed = 0, grid_size = 0;
	long long n_par = 0, n_timesteps = 0;
	particle_t* par_array = NULL;
    Gnode* grid = NULL;
    Vector2D output;
    
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
		fprintf(stderr, "Error in array allocation, simpar.c1!\n");
		exit(1);
	}	
    
    grid = (Gnode*)malloc(grid_size*grid_size*sizeof(Gnode));
	if(grid == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar.c!\n");
		exit(1);
	}	
    
	init_particles(seed, grid_size, n_par, par_array);
    
    for(n = 0; n < n_timesteps; n++)
    {
        //fill grid with zeros
        clearGrid(grid, grid_size);
        //calculate centermass of each cell
        calculate_centermass(grid, grid_size, par_array, n_par);
        compute_newstate(grid, grid_size, par_array, n_par);
    }
    
    output = calculate_finalcentermass(par_array, n_par);
    
    printf("\n%.2lf %.2lf", par_array[0].x, par_array[0].y);
	printf("\n%.2lf %.2lf\n", output.x, output.y);
    
	free(par_array);
    free(grid);
    
	return 0;
}













