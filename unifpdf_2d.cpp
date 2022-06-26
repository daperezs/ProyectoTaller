//
// Created by david on 09/05/2022.
//

#include "unifpdf_2d.h"

//unifpdf_2d
//Last modified 21 November 2013
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

//A 2D uniform PDF function.
//xrange and yrange are inclusive limits of region
//z is a 2D sample point.
//This is used in the clutter function (i.e. we assume a uniform
//distribution of clutter).
//This implementation is fairly standard for a uniform distribution in X and
//Y. This can be problematic for our application if care is not taken.
//With a uniform distribution, we assume there is ZERO probability outside of the region of interest.
//This implies that any measurement outside the region of interest has ZERO probability of
//clutter, when it should have 100% probability of clutter (i.e. we cannot make good measurements
//outside that region).
//For this implementation we never have any measurements outside the region
//of interest, and the state estimates are never predicted/updated outside of that region
//but this might not be the case with other simulations or actual data. BEWARE!

double unifpdf_2d(int xrange[], int yrange[], double z[]){

    int minX, maxX, minY, maxY;
    double val, evalX, evalY;

    minX = xrange[0];
    maxX = xrange[1];
    minY = yrange[0];
    maxY = yrange[1];
    evalX = z[0];
    evalY = z[1];

    if(evalX < minX)
        val = 0;
    else if(evalX > maxX)
        val = 0;
    else if (evalY < minY)
        val = 0;
    else if(evalY > maxY)
        val = 0;
    else
        val = 1 / ((maxX - minX) * (maxY - minY));
    return val;
}
