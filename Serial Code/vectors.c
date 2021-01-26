#include "vectors.h"

Vector2D createVector(double initial_x, double initial_y)
{
    Vector2D v;
    v.x = initial_x;
    v.y = initial_y;
    
    return v;
}

Vector2D addVectors(Vector2D a, Vector2D b)
{
    Vector2D v;
    
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    
    return v;
}

Vector2D subVectors(Vector2D a, Vector2D b)
{
    Vector2D v;
    
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    
    return v;
}

Vector2D normalizeVector(Vector2D v)
{
    Vector2D v2;
    double lenght =sqrt(v.x*v.x + v.y*v.y);
    
    v2.x = v.x/lenght;
    v2.y = v.y/lenght;
    
    return v2;
}

Vector2D scaleVector(Vector2D v, double scale_factor)
{
    Vector2D v2;
    
    v2.x = scale_factor*v.x;
    v2.y = scale_factor*v.y;
    
    return v2;
}

double getLenghtVector(Vector2D v)
{
    return sqrt(v.x*v.x + v.y*v.y);
}

double innerProductVectors(Vector2D a, Vector2D b)
{
    return (a.x*b.x + a.y*b.y);
}












