#ifndef GRID_H
#define GRID_H

#include "particle.h"

#define index(row, column, grid_size) (column + (row*grid_size))

typedef struct _grid_node
{
    double centermass_x, centermass_y;
    double total_mass;
}Gnode;

void clearGrid(Gnode* grid, long grid_size);
void calculate_centermass(Gnode* grid, long grid_size, ParArray* par_array, long long n_par);
Vector2D calculate_force(particle_t par, Gnode* grid, long grid_size);
void compute_newstate(Gnode* grid, long grid_size, ParArray* par_array, long long n_par, ProcParams params, MPI_Comm grid_comm, int *neighbours, int tag, MPI_Datatype mpi_particle, ParArray* send_recvs, ParArray* recv_array, MPI_Comm active_comm);
#endif
