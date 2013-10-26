/* Prototype functions for Mueller Potential */

#include "Structures.h"

#ifndef MuellerPotential_H

#define MuellerPotential_H

void initialiseFourParameter(double* );

void initialiseTwoParameter(double* );

void initialiseA(double* );

void initialiseaa(double* );

void initialisebb(double* );

void initialisecc(double* );

void initialisex0(double* );

void initialisey0(double* );

void initialiser(double* );

void initialisev(double* );

double getMuellerPotential(point );

void getMuellerForce(point* );

#endif