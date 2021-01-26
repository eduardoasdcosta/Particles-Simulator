#include "grid.h"

void clearGrid(Gnode* grid, long grid_size)
{
    int i = 0;
    long size = grid_size*grid_size;
    
    for(i = 0; i < size; i++)
    {
        grid[i].centermass_x = 0.0f;
        grid[i].centermass_y = 0.0f;
        grid[i].total_mass = 0.0f;
    }
}

void calculate_centermass(Gnode* grid, long grid_size, particle_t* par_array, long long n_par)
{
    long long i = 0;
    long j = 0;
    long size = grid_size*grid_size;
    
    for(i = 0; i < n_par; i++)
    {
        j = index(par_array[i].cell_line, par_array[i].cell_column, grid_size);
        grid[j].centermass_x += par_array[i].mass*par_array[i].x;
        grid[j].centermass_y += par_array[i].mass*par_array[i].y;
        grid[j].total_mass += par_array[i].mass;
    }
    
    for(j = 0; j < size; j++)
    {
        if( grid[j].centermass_x != 0)
            grid[j].centermass_x = grid[j].centermass_x / grid[j].total_mass;
        if( grid[j].centermass_y != 0)
            grid[j].centermass_y = grid[j].centermass_y / grid[j].total_mass;
    }
    
}

Vector2D calculate_force(particle_t par, Gnode* grid, long grid_size)
{
    int i = par.cell_line, j = par.cell_column;
    int n1 = 0, n2 = 0;
    int a = 0, b = 0;
    long ix = 0;
    Vector2D resulting_force = createVector(0.0f, 0.0f);
    double ca = 0.0f, co = 0.0f, h = 0.0f, h2 = 0.0f, force = 0.0f, offx = 0.0f, offy = 0.0f;
    
    //calculate resulting force

    for(n1 = -1; n1 < 2; n1++)
    {
        a = i + n1;
        
        if(a < 0)
        {
            a = grid_size - 1;
            offx = -1.0f;
        }
        else if (a >= grid_size)
        {
            a = 0;
            offx = 1.0f;
        }
        else
        {
            offx = 0.0f;
        }
        
        for(n2 = -1; n2 < 2; n2++)
        {
            b = j + n2;
        
            if(b < 0)
            {
                b = grid_size - 1;
                offy = -1.0f;
            }
            else if (a >= grid_size)
            {
                b = 0;
                offy = 1.0f;
            }
            else
            {
                offy = 0.0f;
            }
            
            ix = index(a, b, grid_size);
            //calculate resulting force
            ca = grid[ix].centermass_x + offx - par.x;
            co = grid[ix].centermass_y + offy - par.y;
            h2 = ca*ca + co*co;
            h = sqrt(h2);
            
            if(h < EPSLON)
                continue;
            
            force = G*(par.mass*grid[ix].total_mass)/h2;
            
            resulting_force.x += force*ca/h;
            resulting_force.y += force*co/h;
            
        }
    }
    
    return resulting_force;
    
}

void compute_newstate(Gnode* grid, long grid_size, particle_t* par_array, long long n_par)
{
    Vector2D resulting_force;
    Vector2D acceleration;
    long long i = 0;
    
    for(i = 0; i < n_par; i++)
    {
        resulting_force = calculate_force(par_array[i], grid, grid_size);
        acceleration = scaleVector(resulting_force, (1.0f/par_array[i].mass));
        par_array[i].x += par_array[i].vx + 0.5f*acceleration.x;
        par_array[i].y += par_array[i].vy + 0.5f*acceleration.y;
        par_array[i].vx += acceleration.x;
        par_array[i].vy += acceleration.y;
        
        
        if(par_array[i].x < 0)
            par_array[i].x++;
        else if(par_array[i].x >1)
            par_array[i].x--;
        
        if(par_array[i].y < 0)
            par_array[i].y++;
        else if(par_array[i].y >1)
            par_array[i].y--;
        
       par_array[i].cell_line = (par_array[i].y) * grid_size;
       par_array[i].cell_column = (par_array[i].x) * grid_size;
        
    }
}

Vector2D calculate_finalcentermass(particle_t* par_array, long long n_par)
{
    Vector2D cmass = createVector(0.0f, 0.0f);
    long long i = 0;
    double total_mass = 0.0f;
    
    for(i = 0; i < n_par; i++)
    {
        cmass.x += par_array[i].mass*par_array[i].x;
        cmass.y += par_array[i].mass*par_array[i].y;
        total_mass += par_array[i].mass;
    }
    
    cmass.x /= total_mass;
    cmass.y /= total_mass;
    
    return cmass;
}





















