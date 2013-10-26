/** Utility functions for structures for use with applying Transition Path Sampling to a surface **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Structures.h"
#include "MuellerPotential.h"

void initialisePath(double tStart, double tEnd, double tStep, path* newPath){
   int size;
   
   size = (int)((tEnd - tStart)/tStep) + 1;
   
   (*newPath).points = (point*)malloc(sizeof(point)*(size));
   (*newPath).tStart = tStart;
   (*newPath).tEnd = tEnd;
   (*newPath).tStep = tStep;
}

int pathConcat(path fwdpath, path bwdpath, path* currentPath){
   double t;
   int index;
   /* checks whether tStep is the same for both fwdpath and bwdpath */
   if(bwdpath.tStep != fwdpath.tStep){
      return 0;
   }
   /* checks whether fwdpath begins where bwdpath ends */
   if((bwdpath.points[getPathIndex(bwdpath, bwdpath.tEnd)]).x != (fwdpath.points[getPathIndex(fwdpath, fwdpath.tStart)]).x){
      return 0;
   }
   if((bwdpath.points[getPathIndex(bwdpath, bwdpath.tEnd)]).y != (fwdpath.points[getPathIndex(fwdpath, fwdpath.tStart)]).y){
      return 0;
   }
   
   free((*currentPath).points);
   initialisePath(bwdpath.tStart, fwdpath.tEnd, fwdpath.tStep, currentPath);
   
   for (t = bwdpath.tStart; t <= bwdpath.tEnd; t += bwdpath.tStep){
      index = getPathIndex(bwdpath, t);
      (*currentPath).points[index] = bwdpath.points[index];
   }
   for (t = fwdpath.tStart + fwdpath.tStep; t <= fwdpath.tEnd; t += fwdpath.tStep){
      index = getPathIndex(fwdpath, t);
      (*currentPath).points[index] = fwdpath.points[index];
   }
   return 1;
}

void storeValidPath(path validPath, int storedCount, path** storedPaths){
   path* moreStoredPaths;

   if(storedCount != 1){
      moreStoredPaths = (path*)realloc((*storedPaths), sizeof(path)*storedCount);
	  if (moreStoredPaths != NULL){
		(*storedPaths) = moreStoredPaths;
	  }
   }
   (*storedPaths)[storedCount - 1] = validPath; /* Note at the moment we have no way of making one path equal to another path */
   (*storedPaths)[storedCount - 1] = validPath; /* We will need also to reduce the number of points on the path */
}

int getPathIndex(path currentPath, double time){
   int index;
   
   index = (int)((time - currentPath.tStart)/currentPath.tStep);
   
   return index;
}

point positionOnPath(path currentPath, double time){
   point position;
   int index;
   
   index = getPathIndex(currentPath, time);
   position = currentPath.points[index];
   
   return position;
}

void initialiseSurface(surface* currentSurface){

   (*currentSurface).xIndexMax = (int)(((*currentSurface).width/(*currentSurface).xStep) - 1);
   (*currentSurface).yIndexMax = (int)(((*currentSurface).length/(*currentSurface).yStep) - 1);
   (*currentSurface).height = (double*)malloc(sizeof(double)*(((*currentSurface).xIndexMax + 1) * ((*currentSurface).yIndexMax + 1)));
}

int getSurfaceIndex(surface currentSurface, point currentPoint){
   int index, xIndex, yIndex;
   
   xIndex = (int)((currentPoint.x - currentSurface.startPoint.x)/currentSurface.xStep);
   yIndex = (int)((currentPoint.y - currentSurface.startPoint.y)/currentSurface.yStep);
   index = (yIndex*(currentSurface.xIndexMax + 1)) + xIndex;
   
   return index;
}

void getHeightOnSurface(point* currentPoint, surface* currentSurface){
   int index;
   
   index = getSurfaceIndex(*currentSurface, *currentPoint);
   (*currentSurface).height[index] = getMuellerPotential(*currentPoint);
   (*currentPoint).z = (*currentSurface).height[index];
}