//
// Created by david on 05/05/2022.
//

#include <cstring>
#include "GM_funcionAnonima.h"
#include "unifpdf_2d.h"
/*
//Utility functions
//Used to calculate the indices of two-dimensional target j in a long list of two-dimensional targets
void calculateDataRange2(int j, int out[], int &n)
{
    n=2;
    out[0] = 2*j-1;
    out[1] = 2*j;
}

//Used to calculate the indices of four-dimensional target j in a long list of four-dimensional targets
void calculateDataRange4(int j, int out[], int &n)
{
    n=4;
    out[0] = 4*(j-1) + 1;
    out[1] = 4*(j-1) + 2;
    out[2] = 4*(j-1) + 3;
    out[3] = 4*(j-1) + 4;
}

//Generate clutter function. There are caveats to its use for clutter outside of xrange or yrange - see the comments in unifpdf_2d.m
double clutter_intensity(double z[])
{
    extern int xrange[2], yrange[2];
    extern double lambda_c, V;

    return(lambda_c * V * unifpdf_2d(xrange, yrange, z));
}

void max(double a[][4], int n, double num, double out[][4], int &m)
{
    m=n;
    memcpy(out, a, 4*4*sizeof(double));
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            if(a[i][j] < num) out[i][j] = num;
}
*/

//FALTA: mvnpdf (Matlab)
/*
    birth_intensity = @(x) (0.1 * mvnpdf(x(1:2), birth_mean1(1:2), covariance_birth(1:2,1:2)) + 0.1 * mvnpdf(x(1:2),birth_mean2(1:2), covariance_birth(1:2,1:2)));//Generate birth weight. This only takes into account the position, not the velocity, as Vo&Ma don't say if they use velocity and I assume that they don't. Taken from page 8 of their paper.
    spawn_intensity = @(x, targetState) 0.05 * mvnpdf(x, targetState, covariance_spawn);//Spawn weight, from page 8 of Vo&Ma.
 */
