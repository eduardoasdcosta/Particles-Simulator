#include "grid.h"

#define MASTER_ID 0
#define getId(x,y,size) (y + size*x)


int main(int argc, char** argv)
{
	//simulation related variables
	long seed = 0, grid_size = 0;
	long long n_par = 0, n_timesteps = 0;
	int i = 0, n = 0;
	ParArray par_array, send_arrays[8], recv_array; 
	par_array.index = 0;
	particle_t first_par;
	first_par.array_pos = -1;
	Gnode* grid = NULL;
	double *up = NULL, *down = NULL, *left = NULL, *right = NULL;
	double corners[12];
	int sx = 0, sy = 0, ex = 0, ey = 0;
	double cmx = 0.0f, cmy = 0.0f, cmm = 0.0f;
	double fcmx = 0.0f, fcmy = 0.0f, fcmm = 0.0f;
	//MPI related variables
	int myid, numprocs, namelen, provided, my_grid_id, tag = 0;
	int ndims = 2;
	int dims[2] = {0, 0};
	ProcParams params;
	params.active = 1;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Comm grid_comm, active_comm;
	MPI_Status stat;
	int count = 0, flag = 0;
	int wrap_around[2], coord[2], temp_coord[2], neighbours[8];
	wrap_around[0] = 1; wrap_around[1] = 1;
	//timing variables
	double elapsed_time = 0.0f;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Get_processor_name(processor_name, &namelen);
	MPI_Barrier(MPI_COMM_WORLD);

	if(myid == MASTER_ID)
		elapsed_time = - MPI_Wtime();

	// check if number of arguments is valid
	if(argc != 5)
	{
		fprintf(stderr, "Invalid number of arguments!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	// get arguments from command line
	seed = atol(argv[1]);
	grid_size = atol(argv[2]);
	n_par = atoll(argv[3]);
	n_timesteps = atoll(argv[4]);

	// check if arguments are valid
	if(grid_size < 3 || n_par < 0 || n_timesteps < 1)
	{
		fprintf(stderr, "Invalid arguments!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	
	//create mpi particle struct
	MPI_Datatype mpi_particle;
	int blocklenghts[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype types[9] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_LONG_LONG};
	MPI_Aint offsets[9] = {offsetof(particle_t, x), offsetof(particle_t, y), offsetof(particle_t, vx), offsetof(particle_t, vy), offsetof(particle_t, mass), offsetof(particle_t, cell_line), offsetof(particle_t, cell_column), offsetof(particle_t, active), offsetof(particle_t, array_pos)};
	MPI_Type_create_struct(9, blocklenghts, offsets, types, &mpi_particle);
	MPI_Type_commit(&mpi_particle);
	// create particle array

	if(n_par < 10000000)
	{
		par_array.size = n_par;
		params.array_size = n_par;
	}
	else
	{
		par_array.size = (n_par*1.05f)/numprocs;
		params.array_size = par_array.size*0.1f*8;
	}

	for(i = 0; i < 8; i++)
    {
        send_arrays[i].size = params.array_size;
        send_arrays[i].index = 0;
        send_arrays[i].array = (particle_t*)malloc(send_arrays[i].size*sizeof(particle_t));
        if(send_arrays[i].array == NULL)
        {
            printf("Error in send array malloc\n");
            fflush(stdout);
			MPI_Finalize();
			exit(1);
        }
    }

    recv_array.size = params.array_size;
    recv_array.index = 0;
    recv_array.array = (particle_t*)malloc(recv_array.size*sizeof(particle_t));
    if(recv_array.array == NULL)
    {
        printf("Error in recv array malloc\n");
        fflush(stdout);
		MPI_Finalize();
		exit(1);
    }

	par_array.array = (particle_t*)malloc(par_array.size*sizeof(particle_t));
	if(par_array.array == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	//create grid
	grid = (Gnode*)malloc(grid_size*grid_size*sizeof(Gnode));
	if(grid == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}

	MPI_Dims_create(numprocs, ndims, dims);

	if(dims[0] > dims[1])
	{
		int a = dims[1];
		dims[1] = dims[0];
		dims[0] = a;
	}

	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, 0, &grid_comm);

	MPI_Comm_rank(grid_comm, &my_grid_id);

	MPI_Cart_coords(grid_comm, my_grid_id, 2, coord);

	tag = my_grid_id;
	
	//set number of cells per grid dimension of the subdivision of the grid and save the remainders 
	params.size_x = grid_size/dims[0];
	params.rem_x = grid_size%dims[0];
	params.size_y = grid_size/dims[1];
	params.rem_y = grid_size%dims[1];
	//////////////////////////////////////////////
	if(params.size_x != 0)
	{
		if(coord[0] < params.rem_x)
		{
			params.start_x = coord[0]*params.size_x + coord[0];
			params.end_x = params.start_x + params.size_x + 1;
		}
		else
		{
			params.start_x = coord[0]*params.size_x + params.rem_x;
			params.end_x = params.start_x + params.size_x;
		}
	}
	else
	{
		if(coord[0] < params.rem_x)
		{
			params.start_x = coord[0];
			params.end_x = params.start_x + 1;
			params.size_x = 1;
		}
		else
		{
			params.start_x = 0;
			params.end_x = 0;
		}
	}
	
	if(params.size_y != 0)
	{
		if( coord[1] < params.rem_y)
		{
			params.start_y = coord[1]*params.size_y + coord[1];
			params.end_y = params.start_y + params.size_y + 1;
		}
		else
		{
			params.start_y = coord[1]*params.size_y + params.rem_y;
			params.end_y = params.start_y + params.size_y;
		}
	}
	else
	{
		if(coord[1] < params.rem_y)
		{
			params.start_y = coord[1];
			params.end_y = params.start_y + 1;
			params.size_y = 1;
		}
		else
		{
			params.start_y = 0;
			params.end_y = 0;
		}
	}

	if( (params.start_x == 0 && params.end_x == 0) || (params.start_y == 0 && params.end_y == 0) )
	{
		params.active = 0;
	}
	
	MPI_Comm_split(grid_comm, params.active, my_grid_id, &active_comm);
	/*
	if(myid == MASTER_ID)
		printf("Dims: x %d y %d, Remainders: x %d y %d\n", dims[0], dims[1], params.rem_x, params.rem_y);
	fflush(stdout);
	*/
	//////////////////////////////////////////////
	// initialize particles
	params.grid_id = my_grid_id;
	init_particles(seed, grid_size, n_par, &par_array, params);

	//arrays to get centers of mass from other machines
	up = malloc(3*(params.end_x-params.start_x)*sizeof(double));
	if(up == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	down = malloc(3*(params.end_x-params.start_x)*sizeof(double));
	if(down == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	left = malloc(3*(params.end_y-params.start_y)*sizeof(double));
	if(left == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	right = malloc(3*(params.end_y-params.start_y)*sizeof(double));
	if(right == NULL)
	{
		fprintf(stderr, "Error in array allocation, simpar-mpi.c!\n");
		fflush(stderr);
		MPI_Finalize();
		exit(1);
	}
	
	sx = params.start_x - 1;
	ex = params.end_x + 1;
	sy = params.start_y - 1;
	ey = params.end_y + 1;
	
	if(sx < 0)
		sx = grid_size-1;
	if(ex >= grid_size)
		ex = 0;
	if(sy < 0)
		sy = grid_size-1;
	if(ey >= grid_size)
		ey = 0;
	
	//main loop
	MPI_Barrier(grid_comm); 

	if(params.active)
	{
		//get process neighbour id's and store them
		if(grid_size < dims[0] && coord[0] == 0)
			temp_coord[0] = params.rem_x - 1;
		else
			temp_coord[0] = coord[0] - 1;
		if(grid_size < dims[1] && coord[1] == 0)
			temp_coord[1] = params.rem_y - 1;
		else
			temp_coord[1] = coord[1] - 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[0]);
		temp_coord[0] = coord[0];
		if(grid_size < dims[1] && coord[1] == 0)
			temp_coord[1] = params.rem_y - 1;
		else
			temp_coord[1] = coord[1] - 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[1]);
		if(grid_size < dims[0] && coord[0] == params.rem_x - 1)
			temp_coord[0] = 0;
		else
			temp_coord[0] = coord[0] + 1;
		if(grid_size < dims[1] && coord[1] == 0)
			temp_coord[1] = params.rem_y - 1;
		else
			temp_coord[1] = coord[1] - 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[2]);		
		if(grid_size < dims[0] && coord[0] == 0)
			temp_coord[0] = params.rem_x - 1;
		else
			temp_coord[0] = coord[0] - 1;
		temp_coord[1] = coord[1];
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[3]);
		if(grid_size < dims[0] && coord[0] == params.rem_x - 1)
			temp_coord[0] = 0;
		else
			temp_coord[0] = coord[0] + 1;
		temp_coord[1] = coord[1];
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[4]);
		if(grid_size < dims[0] && coord[0] == 0)
			temp_coord[0] = params.rem_x - 1;
		else
			temp_coord[0] = coord[0] - 1;
		if(grid_size < dims[1] && coord[1] == params.rem_y - 1)
			temp_coord[1] = 0;
		else
			temp_coord[1] = coord[1] + 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[5]);
		temp_coord[0] = coord[0];
		if(grid_size < dims[1] && coord[1] == params.rem_y - 1)
			temp_coord[1] = 0;
		else
			temp_coord[1] = coord[1] + 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[6]);
		if(grid_size < dims[0] && coord[0] == params.rem_x - 1)
			temp_coord[0] = 0;
		else
			temp_coord[0] = coord[0] + 1;
		if(grid_size < dims[1] && coord[1] == params.rem_y - 1)
			temp_coord[1] = 0;
		else
			temp_coord[1] = coord[1] + 1;
		MPI_Cart_rank(grid_comm, temp_coord, &neighbours[7]);
		
		for(n = 0; n < n_timesteps; n++)
		{
			//clear grid and calculate local centers of mass
			clearGrid(grid, grid_size);
			
			calculate_centermass(grid, grid_size, &par_array, n_par);
			
			//share and get centers of mass needed to compute new state
			
			for(i = 0; i < params.end_x-params.start_x; i++)
			{
				up[3*i] = grid[params.start_x + i + (params.start_y)*grid_size].centermass_x;
				up[3*i + 1] = grid[params.start_x + i + (params.start_y)*grid_size].centermass_y;
				up[3*i + 2] = grid[params.start_x + i + (params.start_y)*grid_size].total_mass;
				down[3*i] = grid[params.start_x + i + (params.end_y-1)*grid_size].centermass_x;
				down[3*i + 1] = grid[params.start_x + i + (params.end_y-1)*grid_size].centermass_y;
				down[3*i + 2] = grid[params.start_x + i + (params.end_y-1)*grid_size].total_mass;			
			}
			
			for(i = 0; i < params.end_y-params.start_y; i++)
			{
				left[3*i] = grid[params.start_x + (params.start_y + i)*grid_size].centermass_x;
				left[3*i + 1] = grid[params.start_x + (params.start_y + i)*grid_size].centermass_y;
				left[3*i + 2] = grid[params.start_x + (params.start_y + i)*grid_size].total_mass;
				right[3*i] = grid[params.end_x - 1 + (params.start_y + i)*grid_size].centermass_x;
				right[3*i + 1] = grid[params.end_x - 1 + (params.start_y + i)*grid_size].centermass_y;
				right[3*i + 2] = grid[params.end_x - 1 + (params.start_y + i)*grid_size].total_mass;
			}
			
			MPI_Barrier(active_comm);
			
			MPI_Send(up, 3, MPI_DOUBLE, neighbours[0], tag, grid_comm);
			MPI_Send(up, 3*(params.end_x-params.start_x), MPI_DOUBLE, neighbours[1], tag, grid_comm);
			MPI_Send(&up[3*(params.end_x-params.start_x - 1)], 3, MPI_DOUBLE, neighbours[2], tag, grid_comm);
			MPI_Send(left, 3*(params.end_y-params.start_y), MPI_DOUBLE, neighbours[3], tag, grid_comm);
			MPI_Send(right, 3*(params.end_y-params.start_y), MPI_DOUBLE, neighbours[4], tag, grid_comm);
			MPI_Send(down, 3, MPI_DOUBLE, neighbours[5], tag, grid_comm);
			MPI_Send(down, 3*(params.end_x-params.start_x), MPI_DOUBLE, neighbours[6], tag, grid_comm);
			MPI_Send(&down[3*(params.end_x-params.start_x - 1)], 3, MPI_DOUBLE, neighbours[7], tag, grid_comm);

			MPI_Barrier(active_comm);
			
			MPI_Recv(corners, 3, MPI_DOUBLE, neighbours[0], neighbours[0], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(up, 3*(params.end_x-params.start_x), MPI_DOUBLE, neighbours[1], neighbours[1], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(&corners[3], 3, MPI_DOUBLE, neighbours[2], neighbours[2], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(left, 3*(params.end_y-params.start_y), MPI_DOUBLE, neighbours[3], neighbours[3], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(right, 3*(params.end_y-params.start_y), MPI_DOUBLE, neighbours[4], neighbours[4], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(&corners[6], 3, MPI_DOUBLE, neighbours[5], neighbours[5], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(down, 3*(params.end_x-params.start_x), MPI_DOUBLE, neighbours[6], neighbours[6], grid_comm, MPI_STATUS_IGNORE);
			MPI_Recv(&corners[9], 3, MPI_DOUBLE, neighbours[7], neighbours[7], grid_comm, MPI_STATUS_IGNORE);
	
			grid[sx + sy*grid_size].centermass_x = corners[0];
			grid[sx + sy*grid_size].centermass_y = corners[1];
			grid[sx + sy*grid_size].total_mass = corners[2];
			
			grid[ex + sy*grid_size].centermass_x = corners[3];
			grid[ex + sy*grid_size].centermass_y = corners[4];
			grid[ex + sy*grid_size].total_mass = corners[5];
			
			grid[sx + ey*grid_size].centermass_x = corners[6];
			grid[sx + ey*grid_size].centermass_y = corners[7];
			grid[sx + ey*grid_size].total_mass = corners[8];
			
			grid[ex + ey*grid_size].centermass_x = corners[9];
			grid[ex + ey*grid_size].centermass_y = corners[10];
			grid[ex + ey*grid_size].total_mass = corners[11];
			
			for(i = 0; i < params.end_x-params.start_x; i++)
			{
				grid[params.start_x + i + sy*grid_size].centermass_x = up[3*i] ;
				grid[params.start_x + i + sy*grid_size].centermass_y = up[3*i + 1];
				grid[params.start_x + i + sy*grid_size].total_mass = up[3*i + 2];

				grid[params.start_x + i + ey*grid_size].centermass_x = down[3*i];
				grid[params.start_x + i + ey*grid_size].centermass_y = down[3*i + 1];
				grid[params.start_x + i + ey*grid_size].total_mass = down[3*i + 2];
			}
			
			for(i = 0; i < params.end_y-params.start_y; i++)
			{
				grid[sx + (params.start_y + i)*grid_size].centermass_x = left[3*i] ;
				grid[sx + (params.start_y + i)*grid_size].centermass_y = left[3*i + 1];
				grid[sx + (params.start_y + i)*grid_size].total_mass = left[3*i + 2];
				
				grid[ex + (params.start_y + i)*grid_size].centermass_x = right[3*i];
				grid[ex + (params.start_y + i)*grid_size].centermass_y = right[3*i + 1];
				grid[ex + (params.start_y + i)*grid_size].total_mass = right[3*i + 2];
			}
			//compute new state of particles
			compute_newstate(grid, grid_size, &par_array, n_par, params, grid_comm, neighbours, tag, mpi_particle, send_arrays, &recv_array, active_comm);

		}
		
		MPI_Barrier(active_comm);
		//calculate final centermass
		#pragma omp parallel for private(i) reduction(+:cmx, cmy, cmm) schedule(static) 
		for(i = 0; i < par_array.index; i++)
		{
			if(par_array.array[i].array_pos == 0)
				first_par = par_array.array[i];
			
			cmx += par_array.array[i].x*par_array.array[i].mass;
			cmy += par_array.array[i].y*par_array.array[i].mass;
			cmm += par_array.array[i].mass;
		}
		
		MPI_Reduce(&cmx, &fcmx, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, active_comm);
		MPI_Reduce(&cmy, &fcmy, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, active_comm);
		MPI_Reduce(&cmm, &fcmm, 1, MPI_DOUBLE, MPI_SUM, MASTER_ID, active_comm);
		
		MPI_Barrier(active_comm);
		//send first particle to master
		if(first_par.array_pos != -1 && myid != MASTER_ID)
		{
			MPI_Send(&first_par, 1, mpi_particle, MASTER_ID, myid, MPI_COMM_WORLD);
		}
		MPI_Barrier(active_comm);
		//receive first particle
		if(first_par.array_pos == -1 && myid == MASTER_ID)
		{
			for(i = 0; i < numprocs; i++)
			{
				MPI_Iprobe(i, i, MPI_COMM_WORLD, &flag, &stat);
				if(flag)
				{
					MPI_Get_count(&stat, mpi_particle, &count);
					
					MPI_Recv(&first_par, count, mpi_particle, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
				}
			}
		}

	}
	
	MPI_Barrier(grid_comm); 
	if(myid == MASTER_ID)
	{
		fcmx /= fcmm;
		fcmy /= fcmm;
		printf("%0.2lf %0.2lf\n%0.2lf %0.2lf\n", first_par.x, first_par.y, fcmx, fcmy);
		elapsed_time += MPI_Wtime();
		printf("Process %d of %s: Total execution time is: %0.2lf.\n", my_grid_id, processor_name, elapsed_time);
		fflush(stdout);
	}

	free(up);
	free(down);
	free(left);
	free(right);

	for (int i = 0; i < 8; i++)
    {
        free(send_arrays[i].array);
    }

	free(grid);
	free(par_array.array);
	
	MPI_Finalize();

    return 0;
}
