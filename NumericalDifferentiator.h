/** Protoype functions for Numerical Differentiator **/

#include "Structures.h"

#ifndef NumericalDifferentiator_H

#define NumericalDifferentiator_H

void findlastxPoint(surface , point , point* );

void findlastyPoint(surface , point , point* );

void findnextxPoint(surface , point , point* );

void findnextyPoint(surface , point , point* );

double differentiatex(point , surface );

double differentiatey(point , surface );

#endif