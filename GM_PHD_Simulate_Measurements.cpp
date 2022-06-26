//
// Created by david on 25/05/2022.
//

//GM_PHD_Simulate_Measurements
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

//This file generates simulated measurement data for  the simulation
//described in Vo&Ma.
//There will be gaussian noise on the measurement and Poisson-distributed clutter
//in the environment.

//If you want to use this PHD filter implementation for another problem, you
//will need to replace this script with another one that populates Z,
//zTrue, and simMeasurementHistory (Z is used in a lot of the update code,
//zTrue and simMeasurementHistory are used in GM_PHD_Simulate_Plot)

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <random>
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"
#include "globales.h"

using namespace std;
void matrizVector(double mat[][4], int n, double vec[], int m, double out[], int &nout);
bool isempty(double v[], int n);

void GM_PHD_Simulate_Measurements(void) {
    // Previous define
    extern double dt, simTarget1State[4], simTarget2State[4], F[4][4], simTarget3State[4], k, simTarget3SpawnTime,
            simTarget3Vel[2], xrange[2], yrange[2], prob_detection, sigma_r, noiseScaler, simTarget1History[NUM_F][NUM],
            simTarget2History[NUM_F][NUM], simTarget3History[NUM_F][NUM], clutter[NUM_F][NUM];
    extern int nClutter, c1History, c2History, c3History;

    // New global variables
    extern double clutterX, clutterY, detect1, detect2, detect3, measX1, measY1, measX2, measY2, measX3,
            measY3, randn, Z[2][53], zTrue[2][3];
    int n;
    double out[4];

    std::mt19937 rng;
    std::uniform_real_distribution<float> urd(0, 1);

    //Note: It is possible to get no measurements if the target is not detected
    //and there is no clutter
    cout << "Step Sim: Simulating measurements." << endl;

    //Simulate target movement
    calculate_Jacobian_F(dt, F, n);
    matrizVector(F, 4, simTarget1State, 4, out, n);
    memcpy(simTarget1State, out, 4 * sizeof(double));
    matrizVector(F, 4, simTarget2State, 4, out, n);
    memcpy(simTarget2State, out, 4 * sizeof(double));
    if (!isempty(simTarget3State,4))
    {
        matrizVector(F, 4, simTarget3State, 4, out, n);
        memcpy(simTarget3State, out, 4 * sizeof(double));
    }

    //Spawn target 3 when k = 66
    if(k == simTarget3SpawnTime)
    {
        memcpy(simTarget3State, simTarget1State, 4 * sizeof(double));
        simTarget3State[2] = simTarget3Vel[0];
        simTarget3State[3] = simTarget3Vel[1];
    }

    //Save target movement for plotting
    for(int i = 0; i<4; i++)
    {
        simTarget1History[i][c1History] = simTarget1State[i];
        simTarget2History[i][c2History] = simTarget2State[i];
        simTarget3History[i][c3History] = simTarget3State[i];
    }
    c1History++;
    c2History++;
    c3History++;

    //First, we generate some clutter in the environment.
    memset(clutter,0, NUM_F*NUM*sizeof(double)); //The observations are of the form

    for (int i=0; i<nClutter; i++)
    {
        //variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
        clutterX = urd(rng) * (xrange[1] - xrange[0]) + xrange[0]; //Random number between xrange(1) and xrange(2), uniformly distributed.
        clutterY = urd(rng) * (yrange[1] - yrange[0]) + yrange[0]; //Random number between yrange(1) and yrange(2), uniformly distributed.

        clutter[0][i] = clutterX;
        clutter[1][i] = clutterY;
    }

    //We are not guaranteed to detect the target - there is only a probability

    detect1 = urd(rng);
    detect2 = urd(rng);
    detect3 = urd(rng);

    if(detect1 > prob_detection)
    {
        measX1 = 0;
        measY1 = 0;
    }
    else
    {
        measX1 = simTarget1State[0] + sigma_r * randn * noiseScaler;
        measY1 = simTarget1State[1] + sigma_r * randn * noiseScaler;
    }
    if(detect2 > prob_detection)
    {
        measX2 = 0;
        measY2 = 0;
    }
    else
    {
        measX2 = simTarget2State[0] + sigma_r * randn * noiseScaler;
        measY2 = simTarget2State[1] + sigma_r * randn * noiseScaler;
    }
    if((k >= simTarget3SpawnTime) && (detect3 <= prob_detection))
    {
        measX3 = simTarget3State[0] + sigma_r * randn * noiseScaler;
        measY3 = simTarget3State[1] + sigma_r * randn * noiseScaler;
    }
    else
    {
        measX3 = 0;
        measY3 = 0;
    }


    //Generate true measurement
    memset(Z, 0, 2*53*sizeof(double));
    Z[0][0] = measX1;
    Z[0][1] = measX2;
    Z[0][2] = measX3;
    Z[1][0] = measY1;
    Z[1][1] = measY2;
    Z[1][2] = measY3;

    for(int i = 0; i<2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            zTrue[i][j] = Z[i][j]; //Store for plotting
        }
    }

    //Append clutter
    for(int i = 0; i<2; i++)
    {
        for (int j = 0; j < 50; j++)
        {
            Z[i][j+3] = clutter[i][j]; //Store for plotting
        }
    }

    // FALTA DEFINIR
    //Store history
    //simMeasurementHistory{k} =  Z;

}

void matrizVector(double mat[][4], int n, double vec[], int m, double out[], int &nout)
{
    if(n != m)
    {
        std::cout << "Different dimensions\n";
        exit(0);
    }
    nout = n;

    for(int i=0; i<n; i++)
    {
        out[i] = 0;
        for(int j=0; j<n; j++)
            out[i] += mat[i][j]*vec[j];
    }
}

bool isempty(double v[], int n)
{
    int i = 0;

    while(i<n){
        if(v[i]!=0){
            return false;
        }
        i++;
    }
    return true;
}


