#include "grid.h"

void clearGrid(Gnode* grid, long grid_size)
{
    int i = 0;
    long size = grid_size*grid_size;
    
    // go through every cell of the grid
    for(i = 0; i < size; i++)
    {
        grid[i].centermass_x = 0.0f;
        grid[i].centermass_y = 0.0f;
        grid[i].total_mass = 0.0f;
    }
}

void calculate_centermass(Gnode* grid, long grid_size, ParArray *par_array, long long n_par)
{
    long long i = 0;
    long j = 0;
    long size = grid_size*grid_size;
    
    // go through every particle
    #pragma omp parallel for private(i, j) schedule(static)
    for(i = 0; i < par_array->index; i++)
    {

        // get particle's location on the grid
        j = par_array->array[i].cell_column + par_array->array[i].cell_line*grid_size;

        // compute each particle's contribution to the centermass
        #pragma omp atomic
            grid[j].centermass_x += par_array->array[i].mass*par_array->array[i].x;
        #pragma omp atomic
            grid[j].centermass_y += par_array->array[i].mass*par_array->array[i].y;
        #pragma omp atomic
            grid[j].total_mass += par_array->array[i].mass;
    }
    
    // go through every cell in the grid
    for(j = 0; j < size; j++)
    {
        // compute centermass's position
        if( grid[j].centermass_x != 0)
            grid[j].centermass_x = grid[j].centermass_x / grid[j].total_mass;

        if( grid[j].centermass_y != 0)
            grid[j].centermass_y = grid[j].centermass_y / grid[j].total_mass;
    }
    
}

Vector2D calculate_force(particle_t par, Gnode* grid, long grid_size)
{
    int par_row = par.cell_line, par_col = par.cell_column;
    int adj_row = 0, adj_col = 0;
    int row = 0, column = 0;
    long index = 0;
    Vector2D resulting_force = createVector(0.0f, 0.0f);
    double adjacent = 0.0f, opposite = 0.0f, hypotenuse = 0.0f, force = 0.0f, offx = 0.0f, offy = 0.0f;

    for(adj_row = -1; adj_row < 2; adj_row++)
    {
        row = par_row + adj_row;
        
        if(row < 0)
        {
            row = grid_size - 1;
            offx = -1.0f;
        }
        else if (row >= grid_size)
        {
            row = 0;
            offx = 1.0f;
        }
        else
        {
            offx = 0.0f;
        }
        
        for(adj_col = -1; adj_col < 2; adj_col++)
        {
            column = par_col + adj_col;
        
            if(column < 0)
            {
                column = grid_size - 1;
                offy = -1.0f;
            }
            else if (column >= grid_size)
            {
                column = 0;
                offy = 1.0f;
            }
            else
            {
                offy = 0.0f;
            }
            
            // get particle's location on the grid
            index = column + row*grid_size;

            // compute adjacent and opposite sides
            adjacent = grid[index].centermass_x + offx - par.x;
            opposite = grid[index].centermass_y + offy - par.y;

            // compute hypotenuse (= length)
            hypotenuse = sqrt(adjacent*adjacent + opposite*opposite);
            
            // if distance is below epsilon, proceed to next 
            if(hypotenuse < EPSLON)
                continue;
            
            // compute force from each cell
            force = G*(par.mass*grid[index].total_mass)/(hypotenuse*hypotenuse);

            // add contribution to the resulting force          
            resulting_force.x += force*adjacent/hypotenuse;
            resulting_force.y += force*opposite/hypotenuse;
            
        }
    }
    
    return resulting_force;
    
}

void compute_newstate(Gnode* grid, long grid_size, ParArray* par_array, long long n_par, ProcParams params, MPI_Comm grid_comm, int *neighbours, int tag, MPI_Datatype mpi_particle, ParArray* send_arrays, ParArray* recv_array, MPI_Comm active_comm)
{
    Vector2D resulting_force;
    Vector2D acceleration;
    long long i = 0;
    int j = 0, count = 0;
	int neighbour_id = 0;
    MPI_Status stat;

    for(i = 0; i < 8; i++)
	{
        send_arrays[i].index = 0;
	}
    recv_array->index = 0;

    #pragma omp parallel for private(i, resulting_force, acceleration) schedule(static)
    for(i = 0; i < par_array->index; i++)
    {
        // compute resulting force
        resulting_force = calculate_force(par_array->array[i], grid, grid_size);

        // compute acceleration
        acceleration = scaleVector(resulting_force, (1.0f/par_array->array[i].mass));
  
        // update particle's information
        par_array->array[i].x += par_array->array[i].vx + 0.5f*acceleration.x;
        par_array->array[i].y += par_array->array[i].vy + 0.5f*acceleration.y;
        par_array->array[i].vx += acceleration.x;
        par_array->array[i].vy += acceleration.y;
        
        
        if(par_array->array[i].x < 0)
            par_array->array[i].x++;
        else if(par_array->array[i].x >1)
            par_array->array[i].x--;
        
        if(par_array->array[i].y < 0)
            par_array->array[i].y++;
        else if(par_array->array[i].y >1)
            par_array->array[i].y--;
        
        // update particle's location on the grid
        par_array->array[i].cell_line = (par_array->array[i].y) * grid_size;
        par_array->array[i].cell_column = (par_array->array[i].x) * grid_size;
        
		if(par_array->array[i].cell_line < params.start_y || par_array->array[i].cell_line >= params.end_y || par_array->array[i].cell_column < params.start_x || par_array->array[i].cell_column >= params.end_x)
        {
            par_array->array[i].active = 0;
        }
    }
    
    i = 0;
 
    while(i < par_array->index)
    {
        if(par_array->array[i].active == 0)
        { 
            if(par_array->array[i].cell_column < params.start_x)
            {
                if(par_array->array[i].cell_line < params.start_y)
                {
                    neighbour_id = 0;
                }
                else if(par_array->array[i].cell_line >= params.end_y)
                {
                    neighbour_id = 5;
                }
                else
                {
                    neighbour_id = 3;
                }
            }
            else if(par_array->array[i].cell_column >= params.end_x)
            {
                if(par_array->array[i].cell_line < params.start_y)
                {
                    neighbour_id = 2;
                }
                else if(par_array->array[i].cell_line >= params.end_y)
                {
                    neighbour_id = 7;
                }
                else
                {
                    neighbour_id = 4;
                }
            }
            else
            {
                if(par_array->array[i].cell_line < params.start_y)
                {
                    neighbour_id = 1;
                }
                else if(par_array->array[i].cell_line >= params.end_y)
                {
                    neighbour_id = 6;
                }
            }

            if(send_arrays[neighbour_id].index < send_arrays[neighbour_id].size)
            {
                send_arrays[neighbour_id].array[send_arrays[neighbour_id].index] = par_array->array[i];
                send_arrays[neighbour_id].index = send_arrays[neighbour_id].index + 1;
            }

            par_array->index = par_array->index - 1;
            par_array->array[i] = par_array->array[par_array->index];
            i--;
        }
        
        i++;
        
    }
    
    
    for (int i = 0; i < 8; i++)
    {
        count = 0;
        
		MPI_Barrier(active_comm);
		
        MPI_Send(send_arrays[i].array, send_arrays[i].index, mpi_particle, neighbours[i], tag, grid_comm);
        
        MPI_Probe(neighbours[7-i], neighbours[7-i], grid_comm, &stat);	
		MPI_Get_count(&stat, mpi_particle, &count);

        MPI_Recv(recv_array->array, count, mpi_particle, neighbours[7-i], neighbours[7-i], grid_comm, MPI_STATUS_IGNORE);

        for(j = 0; j < count; j++)
        {
            
            if(recv_array->array[j].array_pos >= n_par)
            {
                printf("got particle with array pos %lld out of bounds\n", recv_array->array[j].array_pos);
                fflush(stdout);
            }
            else
            {
                if(par_array->index < par_array->size)
                {
                    par_array->array[par_array->index] = recv_array->array[j];
                    par_array->index = par_array->index + 1;
                }
            }
            
        }

    }
    /*
    printf("Process %d sent %d particles and received %d particles\n", params.grid_id, tsend, tcount);
    fflush(stdout);
	*/
}





















