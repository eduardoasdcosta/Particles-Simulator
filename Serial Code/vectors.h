#ifndef VECTORS_H
#define VECTORS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _double_vector2D
{
    double x, y;
}Vector2D;


Vector2D createVector(double initial_x, double initial_y);
Vector2D addVectors(Vector2D a, Vector2D b);
Vector2D subVectors(Vector2D a, Vector2D b);
Vector2D normalizeVector(Vector2D v);
Vector2D scaleVector(Vector2D v, double scale_factor);
double getLenghtVector(Vector2D v);
double innerProductVectors(Vector2D a, Vector2D b);

#endif
