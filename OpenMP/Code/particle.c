//.c file for particle related functions
#include "particle.h"

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
{
    long long i = 0;

    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);

        par[i].cell_line = par[i].y * ncside;
        par[i].cell_column = par[i].x * ncside;
    }
}