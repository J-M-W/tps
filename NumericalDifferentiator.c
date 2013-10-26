/** Numerical differentiator for use with Transition Path Sampling on the Mueller surface **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NumericalDifferentiator.h"

void findlastxPoint(surface currentSurface, point currentPoint, point* lastPoint){
   int currentxIndex, currentyIndex, lastxIndex, lastyIndex;
   
   currentxIndex = (int)((currentPoint.x - currentSurface.startPoint.x)/currentSurface.xStep);
   currentyIndex = (int)((currentPoint.y - currentSurface.startPoint.y)/currentSurface.yStep);

   if(currentxIndex == 0){
      lastxIndex = currentxIndex;
      lastyIndex = currentyIndex;
   } else {
      lastxIndex = currentxIndex - 1;
      lastyIndex = currentyIndex;
   }
   
   /* Find the point belonging to the index */
   (*lastPoint).x = (double)((lastxIndex * currentSurface.xStep) + currentSurface.startPoint.x);
   (*lastPoint).y = (double)((lastyIndex * currentSurface.yStep) + currentSurface.startPoint.y);
}

void findlastyPoint(surface currentSurface, point currentPoint, point* lastPoint){
   int currentxIndex, currentyIndex, lastxIndex, lastyIndex;
   
   currentxIndex = (int)((currentPoint.x - currentSurface.startPoint.x)/currentSurface.xStep);
   currentyIndex = (int)((currentPoint.y - currentSurface.startPoint.y)/currentSurface.yStep);

   if(currentyIndex == 0){
      lastxIndex = currentxIndex;
      lastyIndex = currentyIndex;
   } else {
      lastxIndex = currentxIndex;
      lastyIndex = currentyIndex - 1;
   }
   
   /* Find the point belonging to the index */
   (*lastPoint).x = (double)((lastxIndex * currentSurface.xStep) + currentSurface.startPoint.x);
   (*lastPoint).y = (double)((lastyIndex * currentSurface.yStep) + currentSurface.startPoint.y);
}

void findnextxPoint(surface currentSurface, point currentPoint, point* nextPoint){
   int currentxIndex, currentyIndex, nextxIndex, nextyIndex;
   
   currentxIndex = (int)((currentPoint.x - currentSurface.startPoint.x)/currentSurface.xStep);
   currentyIndex = (int)((currentPoint.y - currentSurface.startPoint.y)/currentSurface.yStep);
   
   if(currentxIndex == currentSurface.xIndexMax){
      nextxIndex = currentxIndex;
      nextyIndex = currentyIndex;
   } else {
      nextxIndex = currentxIndex + 1;
      nextyIndex = currentyIndex;
   }
   
   /* Find the point belonging to the index */
   (*nextPoint).x = (double)((nextxIndex * currentSurface.xStep) + currentSurface.startPoint.x);
   (*nextPoint).y = (double)((nextyIndex * currentSurface.yStep) + currentSurface.startPoint.y);
}

void findnextyPoint(surface currentSurface, point currentPoint, point* nextPoint){
   int currentxIndex, currentyIndex, nextxIndex, nextyIndex;
   
   currentxIndex = (int)((currentPoint.x - currentSurface.startPoint.x)/currentSurface.xStep);
   currentyIndex = (int)((currentPoint.y - currentSurface.startPoint.y)/currentSurface.yStep);
   
   if(currentyIndex == currentSurface.yIndexMax){
      nextxIndex = currentxIndex;
      nextyIndex = currentyIndex;
   } else {
      nextxIndex = currentxIndex;
      nextyIndex = currentyIndex + 1;
   }
   
   /* Find the point belonging to the index */
   (*nextPoint).x = (double)((nextxIndex * currentSurface.xStep) + currentSurface.startPoint.x);
   (*nextPoint).y = (double)((nextyIndex * currentSurface.yStep) + currentSurface.startPoint.y);
}

double differentiatex(point currentPoint, surface currentSurface){   
   point lastPoint, nextPoint;
   double gradient;
   
   /* Find the last and nextPoint */
   findlastxPoint(currentSurface, currentPoint, &lastPoint);
   findnextxPoint(currentSurface, currentPoint, &nextPoint);
    
   /* Get the z co-ordinates of next and lastPoint */
   getHeightOnSurface(&lastPoint, &currentSurface);
   getHeightOnSurface(&nextPoint, &currentSurface);
   
   /* Evaluate the gradient of a line drawn between next and lastPoint at constant y */
   gradient = (nextPoint.z - lastPoint.z)/(nextPoint.x - lastPoint.x);
   
   return gradient;
}

double differentiatey(point currentPoint, surface currentSurface){
   point lastPoint, nextPoint;
   double gradient;
   
   /* Find the last and nextPoint */
   findlastyPoint(currentSurface, currentPoint, &lastPoint);
   findnextyPoint(currentSurface, currentPoint, &nextPoint);
   
   /* Get the z co-ordinates of next and lastPoint */
   getHeightOnSurface(&lastPoint, &currentSurface);
   getHeightOnSurface(&nextPoint, &currentSurface);
   
   /* Evaluate the gradient of a line drawn between next and lastPoint at constant x */
   gradient = (nextPoint.z - lastPoint.z)/(nextPoint.y - lastPoint.y);
   
   return gradient;
}