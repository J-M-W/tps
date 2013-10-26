#include "structures.h"
#include <stdlib.h>

int index2D(int indX, int indY, int nX)
{
    /* Converts 2D indices (indX, indY) into the equivalent index for the flattened array.
    2D array must be of x-length nX */
    return indX + nX*indY;
}

double interpolate(surface* surf, double* r, int order)
{
    /* Interpolates the values on 'surf' at the point 'r'
    order=0 - nearest-neighbour
    order=1 - bilinear
    order=2 - bicubic */

    double energy=0;
    double ind[2]; /* surface-style indices of point 'r' (not constrained to integer values) */

    ind[0] = (r[0] - (surf->startPoint).x)/surf->xStep;
    ind[1] = (r[1] - (surf->startPoint).y)/surf->yStep;
    if(order==1)
    {
        double tl, tr, bl, br; /* Top-left, top-right, bottom-left, bottom-right.
                                  Bottom and right are in the direction of increasing x and y */
        double x1, x2, y1, y2;

        /* Calculate (x, y) coordinates of grid nodes */
        x1 = surf->startPoint.x + (int)ind[0]*surf->xStep;
        x2 = surf->startPoint.x + ((int)ind[0]+1)*surf->xStep;
        y1 = surf->startPoint.y + (int)ind[0]*surf->yStep;
        y2 = surf->startPoint.y + ((int)ind[0]+1)*surf->yStep;

        tl = surf->height[index2D((int)ind[0], (int)ind[1]+1, surf->xIndexMax+1)];
        tr = surf->height[index2D((int)ind[0]+1, (int)ind[1]+1, surf->xIndexMax+1)];
        bl = surf->height[index2D((int)ind[0], (int)ind[1], surf->xIndexMax+1)];
        br = surf->height[index2D((int)ind[0]+1, (int)ind[1], surf->xIndexMax+1)];

        energy=0;
        energy += bl*(x2-r[0])*(y2-r[1]);
        energy += br*(r[0]-x1)*(y2-r[1]);
        energy += tl*(x2-r[0])*(r[1]-y1);
        energy += tr*(r[0]-x1)*(r[1]-y1);
        energy /= (x2-x1)*(y2-y1);
        return energy;
    }
    else if(order==2)
    {
        puts("Bicubic interpolation not yet implemented.");
        exit(EXIT_FAILURE);
    }
    else if(order==0)
    {
        puts("Nearest-neighbour interpolation not yet implemented");
    }
    else
    {
        puts("Order most be 0, 1 or 2");
    }

    /* Get neighbours */

}
