/** Framework structure of code to perform Transition Path Sampling on the Mueller surface **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Structures.h"
#include "NumericalDifferentiator.h"
#include "mt19937ar.h"

#define TRUE 1
#define FALSE 0

void generateEnergySurface(surface* energySurface){
   /* set values for the width and length of the surface */
   /* choose xStep and yStep */
   initialiseSurface(energySurface);
}

void findStationaryPoints(surface energySurface, point* stationaryPoints){
   point currentPoint, lastPoint, nextPoint;
   point* moreStationaryPoints;
   int yindex, xindex, index;
   double xgradient, ygradient;
   int y, x, SPcount;
   
   for(y = 0; y <= energySurface.yIndexMax; y++){
      yindex = y;
	  for(x = 0; x <= energySurface.xIndexMax; x++){
		  xindex = x;

		  currentPoint.x = (double)((xindex * energySurface.xStep) + energySurface.startPoint.x);
		  currentPoint.y = (double)((yindex * energySurface.yStep) + energySurface.startPoint.y);

		  xgradient = differentiatex(currentPoint, energySurface);

		  ygradient = differentiatey(currentPoint, energySurface);

		  if((xgradient == 0) && (ygradient == 0)){
			  SPcount += 1;
			  if(SPcount != 1){
				  moreStationaryPoints = (point*)realloc(stationaryPoints, sizeof(point)*SPcount);
				  if (moreStationaryPoints != NULL){
					stationaryPoints = moreStationaryPoints;
				  }
			  }
			  stationaryPoints[SPcount - 1].x = currentPoint.x;
			  stationaryPoints[SPcount - 1].y = currentPoint.y;
		  }
	  }
   }
}

void findMinima(){
	/* Stuff with the Hessian - fun times await! */
   
}

void setEndRegion(point minimum, path* regionPerimeter){
   /* Define a circle around each analytically determined minima to be the region perimeter */
   initialisePath(tStart, tEnd, tStep, regionPerimeter);
}

int transitionPathSampling(){
   /* apply Transition Path Sampling */
   path currentPath;
   path* storedPaths;
   point* stationaryPoints;
   surface energySurface;
   point startPoint, endPoint;
   int storedCount, storedMax /* Some predefined maximum number of paths to be stored */;

   stationaryPoints = (point*)malloc(sizeof(point)*1);
   storedPaths = (path*)malloc(sizeof(path)*1);

   findStationaryPoints(energySurface, stationaryPoints);
   
   /* Decide values for tStart, tEnd and tStep */
   currentPath.tStart = 0; /* tStart set to 0 for now */
   currentPath.tEnd = 10; /* tEnd set to 10 for now */ /* These values need to be determined! */
   currentPath.tStep = 1; /* tStep set to 1 for now */

   if(getInitialPath(energySurface, &currentPath) == 0){
	   return EXIT_FAILURE;
   }
   
   storedCount = 0;

   while (storedCount < storedMax){
      trialMove(&energySurface, &currentPath, storedCount, &storedPaths);
   }


   free(stationaryPoints);
   free(currentPath.points);
   free(energySurface.height);
      
}

int getInitialPath(surface energySurface, path* currentPath){
   point minimumA, minimumB;
   double t;
   int index, nSteps /* the number of time steps over the whole path */, xStep, yStep;

   initialisePath((*currentPath).tStart, (*currentPath).tEnd, (*currentPath).tStep, currentPath);

   /* Calls find minima and so knows the x and y co-ordinates of the desired start and end points */

   (*currentPath).startPoint.x = minimumA.x;
   (*currentPath).startPoint.y = minimumA.y;

   nSteps = ((*currentPath).tEnd - (*currentPath).tStart)/(*currentPath).tStep;
   xStep = ((*currentPath).endPoint.x - (*currentPath).startPoint.x)/nSteps;
   yStep = ((*currentPath).endPoint.y - (*currentPath).startPoint.y)/nSteps;

   for(t = (*currentPath).tStart; t < (*currentPath).tEnd; t += (*currentPath).tStep){
      index = getPathIndex((*currentPath), t);
	  (*currentPath)[index].x = (*currentPath).tStart.x + t*xStep;
	  (*currentPath)[index].y = (*currentPath).tStart.y + t*yStep;
	  (*currentPath)[index].z = getHeightOnSurface((*currentPath)[index], &energySurface);
   }

   /* Validate that the path ends in the correct place */

   if((*currentPath)[getPathIndex((*currentPath), (*currentPath).tEnd)].x != minimumB.x){
      return 0;
   } else if ((*currentPath)[getPathIndex((*currentPath), (*currentPath).tEnd)].y != minimumB.y){
	  return 0;
   } else {
	  return 1;
   }
}

void trialMove(surface* energySurface, path* currentPath, int storedCount, path** storedPaths){
   path fwdpath, bwdpath;
   point slicePoint;
   double sliceTime;
   int trialMoveValid;

   while(trialMoveValid == TRUE){
      getRandomSlice(*currentPath, &sliceTime, &slicePoint);
      if(addDisplacement((*energySurface), &slicePoint) == FALSE) {
         trialMoveValid = FALSE;
         break;
      }
      if(integrateForward((*energySurface), slicePoint, sliceTime, &fwdpath) == FALSE) {
         trialMoveValid = FALSE;
         break;
      }
      if (integrateBackward((*energySurface), slicePoint, sliceTime, &bwdpath) == FALSE) {
         trialMoveValid = FALSE;
         break;
      }
      
      pathConcat(fwdpath, bwdpath, currentPath);
	  storedCount++;
      storeValidPath(*currentPath, storedCount, storedPaths);
   }
}

void getRandomSlice(path currentPath, double* sliceTime, point* slicePoint){
   /* generate a random number within the bounds of the path's lifetime */
   /* function 'genrand_real1()' is from mt19937ar.c; generates a random number between 0 and 1 inclusive */
   
   *sliceTime = (genrand_real1() * (currentPath.tEnd - currentPath.tStart)) + currentPath.tStart;
   
   /* use this to choose a random slice of the trajectory */
   
   *slicePoint = positionOnPath(currentPath, *sliceTime);
}

int addDisplacement(surface energySurface, point* slicePoint){
   double displacementx, displacementy;
   /* generate a random number and assign this to displacement */
   
   displacementx = genrand_real1() * (energySurface.width * 0.001); /* Random number is chosen to be arbitrarily small */
   displacementy = genrand_real1() * (energySurface.length * 0.001);
   
   /* modify the time slice */
   
   (*slicePoint).x += displacementx;
   (*slicePoint).y += displacementy;
   
   if(isPointInSurfaceBounds((*slicePoint), energySurface)){
      /* assign (*slicePoint).z to the height on the surface at slicePoint.x, slicePoint.y */
      getHeightOnSurface(slicePoint, &energySurface);
      return TRUE;
   } else {
      return FALSE;
   }
   
}

int isPointInSurfaceBounds(point slicePoint, surface energySurface){
   if((slicePoint.x < 0) || (slicePoint.x >= energySurface.width)){
      return 0;
   } if ((slicePoint.y < 0) || (slicePoint.y >= energySurface.length)){
      return 0;
   } else {
      return 1;
   }
}

int integrateForward(surface energySurface, point slicePoint, double sliceTime, path regionBPerimeter, path* fwdpath){
   /* integrate time forward from the random time slice */
   
   /* Verlet algorithm */
   
   point endPoint;
   
   endPoint.x = ((*fwdpath).points[(*fwdpath).tEnd]).x;
   endPoint.y = ((*fwdpath).points[(*fwdpath).tEnd]).y;
   
   return validEndPoint(endPoint, energySurface, (*fwdpath).tEnd, regionBPerimeter);
}

int validEndPoint(point endPoint, surface currentSurface, double endTime, path regionPerimeter){
   /* checks if the end point of the path will lie in the desired end region */
   int intersections, xIndex;
   
   xIndex = (int)(endPoint.x/currentSurface.xStep);
   currentSurface.xIndexMax;
   
   /** 
   
   for(i = xIndex; i < currentSurface.xIndexMax; i++){
      if( point lies on the perimeter ){
         for(j = yIndex; j < currentSurface.yIndexMax; j++){
            if( point lies on the perimeter ){
               intersections++;
            }
         }
      }
   }
   
   **/
   
   if(intersections % 2){
      return 0;
   } else {
      return 1;
   }
}

int integrateBackward(surface energySurface, point slicePoint, double sliceTime, path regionAPerimeter, path* bwdpath){
   /* integrate time backwards from the slice time */
   
   /* Verlet algorithm */
   
   point endPoint;
   
   endPoint.x = ((*bwdpath).points[(*bwdpath).tStart]).x;
   endPoint.y = ((*bwdpath).points[(*bwdpath).tStart]).y;
   
   return validEndPoint(endPoint, energySurface, (*bwdpath).tStart, regionAPerimeter);
}
