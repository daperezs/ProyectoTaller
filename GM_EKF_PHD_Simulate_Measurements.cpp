//
// Created by david on 02/06/2022.
//

/*
%GM_PHD_Simulate_Measurements
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

%This file generates simulated measurement data for the linear KF simulation
%described in Vo&Ma. There is a nonlinear simulation described in the paper
%but I have not gotten around to implementing it.
%The sensor measurements are range-bearing measurements, similar to those
%from an FMCW radar. The expected error magnitude is independent of
%measurement (i.e. range error is not proportional to range, bearing error
%is not proportional to bearing).

%There will be gaussian noise on the measurement and Poisson-distributed clutter
%in the environment.
%Note: It is possible to get no measurements if the target is not detected
%and there is no clutter

%If you want to use this PHD filter implementation for another problem, you
%will need to replace this script with another one that populates Z,
%zTrue, and simMeasurementHistory (Z is used in a lot of the update code,
%zTrue and simMeasurementHistory are used in GM_PHD_Simulate_Plot, and
%simMeasurementHistory is used in GM_EKF_PHD_Create_Birth)
 */

#include <iostream>
#include <random>
#include <cstring>
#include "GM_EKF_PHD_Simulate_Measurements.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"

using namespace std;

void matrizVectorEKF(double mat[][4], int n, double vec[], int m, double out[], int &nout);
bool isemptyEKF(double v[], int n);
void hEKF(double xS, double yS, double hS, double xL, double yL, double out[], int &n);

void GM_EKF_PHD_Simulate_Measurements(void)
{
    //Previous define variables
    extern double F[4][4], dt, simTarget1State[4], simTarget2State[4], simTarget3State[4], k, simTarget3SpawnTime,
            simTarget3Vel[2], simTarget1History[NUM_F][NUM], simTarget2History[NUM_F][NUM], simTarget3History[NUM_F][NUM],
            xrange[2], yrange[2], clutterX, clutterY, x_sensor[3], detect1, detect2, detect3, prob_detection, sigma_r,
            noiseScaler, randn, sigma_theta, Z[2][53], zTrue[2][3];;
    extern int c1History, c2History, c3History, nClutter;

    //New variables
    extern double clutterRTheta[2], measR1, measTheta1, simMeas1[2], measR2, measTheta2, simMeas2[2], measR3, measTheta3,
            simMeas3[2], clutter[NUM_F][NUM];

    int n, m;
    double out[4];

    std::mt19937 rng;
    std::uniform_real_distribution<float> urd(0, 1);

    //There will be gaussian noise on the measurement and Poisson clutter
    //in the environment.
    //Note: It is possible to get no measurements if the target is not detected
    //and there is no clutter
    cout << "Step Sim: Simulating measurements." << endl;

    //Simulate target movement
    //F = calculate_Jacobian_F(dt);
    calculate_Jacobian_F(dt, F, n);
    matrizVectorEKF(F, 4, simTarget1State, 4, out, n);
    memcpy(simTarget1State, out, 4 * sizeof(double));
    matrizVectorEKF(F, 4, simTarget2State, 4, out, n);
    memcpy(simTarget2State, out, 4 * sizeof(double));
    if (!isemptyEKF(simTarget3State,4))
    {
        matrizVectorEKF(F, 4, simTarget3State, 4, out, n);
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

        h(x_sensor[0], x_sensor[1], x_sensor[2], clutterX, clutterY,clutterRTheta,m);

        clutter[0][i] = clutterRTheta[0];
        clutter[1][i] = clutterRTheta[1];
    }

    //We are not guaranteed to detect the target - there is only a probability
    //that we will, controlled by prob_detection
    detect1 = urd(rng);
    detect2 = urd(rng);
    detect3 = urd(rng);
    //Target 1
    if(detect1 > prob_detection)
    {
        measR1 = 0;
        measTheta1 = 0;
    }else
    {
        h(x_sensor[0],x_sensor[1], x_sensor[2], simTarget1State[0], simTarget1State[1], simMeas1, m); //Generate measurement
        measR1 = simMeas1[0] + sigma_r * randn * noiseScaler; //Add gaussian noise. We could add the noise at the measurement-generating stage, by adding it to the simulated target state.
        measTheta1 = simMeas1[1] + sigma_theta * randn * noiseScaler; //Add gaussian noise. We could add the noise at the measurement-generating stage, by adding it to the simulated target state.
    }
    //Target 2
    if(detect2 > prob_detection)
    {
        measR2 = 0;
        measTheta2 = 0;
    }else
    {
        h(x_sensor[0], x_sensor[1], x_sensor[2],  simTarget2State[0], simTarget2State[1], simMeas2, m);
        measR2 = simMeas2[0] + sigma_r * randn * noiseScaler;
        measTheta2 = simMeas2[1] + sigma_theta * randn * noiseScaler;
    }
    //Target 3
    if((k >= simTarget3SpawnTime) && (detect3 <= prob_detection))
    {
        h(x_sensor[0], x_sensor[1], x_sensor[2],  simTarget3State[0], simTarget3State[1], simMeas3, m);
        measR3 = simMeas3[0] + sigma_r * randn * noiseScaler;
        measTheta3 = simMeas3[1] + sigma_theta * randn * noiseScaler;
    }else
    {
        measR3 = 0;
        measTheta3 = 0;
    }

    //Generate true measurement
    memset(Z, 0, 2*53*sizeof(double));
    Z[0][0] = measR1;
    Z[0][1] = measR2;
    Z[0][2] = measR3;
    Z[1][0] = measTheta1;
    Z[1][1] = measTheta2;
    Z[1][2] = measTheta3;

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
    //FALTA
    //simMeasurementHistory{k} =  Z;
}

void matrizVectorEKF(double mat[][4], int n, double vec[], int m, double out[], int &nout)
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

bool isemptyEKF(double v[], int n)
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

//h = @(xS,yS,hS,xL,yL)[hypot(xL-xS,yL-yS);atan2(yL-yS,xL-xS)-hS]
void hEKF(double xS, double yS, double hS, double xL, double yL, double out[], int &n)
{
    n=2;
    out[0]=hypot(xL-xS,yL-yS);
    out[1]=atan2(yL-yS,xL-xS)-hS;
}

