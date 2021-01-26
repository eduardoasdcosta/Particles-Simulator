//.c file for particle related functions
#include "particle.h"

void init_particles(long seed, long ncside, long long n_part, ParArray *par, ProcParams params)
{
    long long i = 0;
    particle_t temp;
    srandom(seed);
	
    for(i = 0; i < n_part; i++)
    {
        temp.x = RND0_1;
        temp.y = RND0_1;
        temp.vx = RND0_1 / ncside / 10.0;
        temp.vy = RND0_1 / ncside / 10.0;
        temp.mass = RND0_1 * ncside / (G * 1e6 * n_part);

        temp.cell_line = temp.y * ncside;
        temp.cell_column = temp.x * ncside;

        if( (temp.cell_column < params.start_x || temp.cell_column >= params.end_x) || (temp.cell_line < params.start_y || temp.cell_line >= params.end_y) )
        {
        }
        else
        {
            temp.active = 1;
            temp.array_pos = i;
            par->array[par->index] = temp;
            
            par->index = par->index + 1;
        }
    }
}
