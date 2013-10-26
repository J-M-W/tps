/** Routines for computing the potential and force due to the Mueller Potential **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MuellerPotential.h"

void initialiseFourParameter(double* fourParameter){
   /* initialises the Mueller parameters that have 4 elements in the array */
   *(fourParameter) = (double*)malloc(sizeof(double)*(4));
}

void initialiseTwoParameter(double* twoParameter){
   /* initialises the Mueller parameters that have 2 elements in the array */
   *(twoParameter) = (double*)malloc(sizeof(double)*(4));
}

void initialiseA(double* MuellerA){
   initialiseFourParameter(MuellerA);
   /* Fill the array with its predefined values */
   (*MuellerA)[0] = -200.0;
   (*MuellerA)[1] = -100.0;
   (*MuellerA)[2] = -170.0;
   (*MuellerA)[3] = 15.0;
}

void initialiseaa(double* Muelleraa){
   initialiseFourParameter(Muelleraa);
   (*Muelleraa)[0] = -1.0;
   (*Muelleraa)[1] = -1.0;
   (*Muelleraa)[2] = -6.5;
   (*Muelleraa)[3] = 0.7;
}

void initialisebb(double* Muellerbb){
   initialiseFourParameter(Muellerbb);
   (*Muellerbb)[0] = 0;
   (*Muellerbb)[1] = 0;
   (*Muellerbb)[2] = 11.0;
   (*Muellerbb)[3] = 0.6;
}

void initialisecc(double* Muellercc){
   initialiseFourParameter(Muellercc);
   (*Muellercc)[0] = -10.0;
   (*Muellercc)[1] = -10.0;
   (*Muellercc)[2] = -6.5;
   (*Muellercc)[3] = 0.7;
}

void initialisex0(double* Muellerx0){
   initialiseFourParameter(Muellerx0);
   (*Muellerx0)[0] = 1.0;
   (*Muellerx0)[1] = 0.0;
   (*Muellerx0)[2] = -0.5;
   (*Muellerx0)[3] = -1.0;
}

void initialisey0(double* Muellery0){
   initialiseFourParameter(Muellery0);
   (*Muellery0)[0] = 0.0;
   (*Muellery0)[1] = 0.5;
   (*Muellery0)[2] = 1.5;
   (*Muellery0)[3] = 1.0;
}

void initialiser(double* Muellerr){
   initialiseTwoParameter(Muellerr);
   (*Muellerr)[0] = 0.5;
   (*Muellerr)[1] = 0.0;
}

void initialisev(double* Muellerv){
   initialiseTwoParameter(Muellerv);
   (*Muellerv)[0] = 0.0;
   (*Muellerv)[1] = 0.0;
}
  
double getMuellerPotential(point currentPoint){
   double potential = 0.0, argument;
   double* MuellerA, Muelleraa, Muellerbb, Muellercc, Muellerx0, Muellery0, Muellerr, Muellerv;
   
   initialiseA(MuellerA);
   initialiseaa(Muelleraa);
   initialisebb(Muellerbb);
   initialisecc(Muellercc);
   initialisex0(Muellerx0);
   initialisey0(Muellery0);
   initialiser(Muellerr);
   initialisev(Muellerv);
   
   int i;
   
   for(i = 0, i < 3, i++){
      argument = Muelleraa[i]*((currentPoint.x - Muellerx0[i])*(currentPoint.x - Muellerx0[i]));
      argument += Muellerbb[i]*(currentPoint.x - Muellerx0[i])*(currentPoint.y - Muellery0[i]);
      argument += Muellercc[i]*((currentPoint.y - Muellery0[i])*(currentPoint.y - Muellery0[i]));
      potential += MuellerA[i]*exp(argument);
   }
   
   return potential;
}

void getMuellerForce(point* currentPoint){
   double force = 0.0, argument;
   double* MuellerA, Muelleraa, Muellerbb, Muellercc, Muellerx0, Muellery0, Muellerr, Muellerv;
   
   initialiseA(MuellerA);
   initialiseaa(Muelleraa);
   initialisebb(Muellerbb);
   initialisecc(Muellercc);
   initialisex0(Muellerx0);
   initialisey0(Muellery0);
   initialiser(Muellerr);
   initialisev(Muellerv);
   
   int i;
   double* xTemp, yTemp;
   
   (*currentPoint).xForce = 0.0;
   (*currentPoint).yForce = 0.0;
   
   for(i = 0, i < 3, i++){
      argument = Muelleraa[i]*((currentPoint.x - Muellerx0[i])*(currentPoint.x - Muellerx0[i]));
      argument += Muellerbb[i]*(currentPoint.x - Muellerx0[i])*(currentPoint.y - Muellery0[i]);
      argument += Muellercc[i]*((currentPoint.y - Muellery0[i])*(currentPoint.y - Muellery0[i]));
      potential += MuellerA[i]*exp(argument);
      
      xTemp = 2*(Muelleraa[i]*(currentPoint.x - Muellerx0[i]));
      xTemp += Muellerbb[i]*(currentPoint.y - Muellery0[i]);
      
      yTemp = 2*(Muellercc[i]*(currentPoint.y - Muellery0[i]);
      yTemp += Muellerbb[i]*(currentPoint.x - Muellerx0[i]);
      
      (*currentPoint).xForce += xTemp*MuellerA[i]*exp(argument);
      (*currentPoint).yForce += yTemp*MuellerA[i]*exp(argument);
      
   }
   
   (*currentPoint).xForce = -(*currentPoint).xForce;
   (*currentPoint).yForce = -(*currentPoint).yForce;
}