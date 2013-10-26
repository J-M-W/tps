/** Structures for use with applying Transition Path Sampling to a surface and the prototypes of their utility functions **/

#ifndef Structures_H

#define Structures_H

typedef struct point{
   double x, y, z;
   double xVelocity, yVelocity;
   double xForce, yForce;
} point;

typedef struct path{
   point* points;
   point startPoint, endPoint;
   double tStart, tEnd, tStep;
} path;

typedef struct surface{
   double* height;
   point startPoint;
   double width, length, xStep, yStep;
   int xIndexMax, yIndexMax;
} surface;

void initialisePath(double , double , double , path* );

int pathConcat(path , path , path* );

void storeValidPath(path , int , path** );

int getPathIndex(path , double );

point positionOnPath(path , double );

void initialiseSurface(surface* );

int getSurfaceIndex(surface , point );

void getHeightOnSurface(point* , surface* );

#endif