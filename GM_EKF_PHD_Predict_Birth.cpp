//
// Created by david
//

/*
%GM_EKF_PHD_Predict_Birth
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

%This file performs prediction for newly birthed and spawned targets
%This is necessary since their positions were initialised in a previous timestep and they may
%have moved between then and now
%We iterate through j = 1..J_b,k where J_b,k is number of birthed landmarks
%at time k (but the landmarks correspond to the measurement from time k-1)
%and j = 1...J_s,k where J_s,k is the number of spawned landmarks at time k
%(but the landmarks correspond to the measurement from time k-1).
%For this implementation, birthed and spawned targets are identical except
%they have weights that are calculated from different functions, and
%different starting covariances. A target is spawned or birthed depending
%on which weighting function will give it a higher starting weight.
*/

#include <iostream>
#include "GM_EKF_PHD_Predict_Birth.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"
#include "GM_PHD_Construct_Update_Components.h"

void matrizxMatrizPB(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc);

using namespace std;


void GM_EKF_PHD_Predict_Birth(void)
{
    //Previos define variables
    extern double m_birth_before_prediction[NUM_F][NUM], m_birth[NUM_F][NUM], F[4][4], Q[4][4], dt, P_birth[NUM_F][NUM],
            P_spawn[NUM_F][NUM], m_spawn[NUM_F][NUM], thisM[4],  w_birth[NUM], w_spawn[NUM];
    extern int numBirthedTargets, numSpawnedTargets, i, VERBOSE, P_range[4];

    //New variables
    extern double G[4][4];

    int n, m;
    double temp[4][4], out[4], temp1[4], out1[4];

    //The means of these will be the position in Cartesian space where they are detected
    //The covariances & weights are calculated according to Vo&Ma
    cout << "Step 1: Prediction for birthed and spawned targets." << endl;

    //Need to store these BEFORE prediction for use in the update step.
    for(int i = 0; i < numBirthedTargets; i++)
    {
        m_birth_before_prediction[0][i] = m_birth[0][i];
        m_birth_before_prediction[1][i] = m_birth[1][i];
        m_birth_before_prediction[2][i] = m_birth[2][i];
        m_birth_before_prediction[3][i] = m_birth[3][i];
    }

    //Perform prediction for birthed targets using birthed velocities.
    for(int j = 0; j <numBirthedTargets; j++)
    {
        i = i + 1;

        calculate_Jacobian_F(dt, F, n);
        calculate_Jacobian_G(dt, G, n);
        //m_birth(:,j) = F * m_birth(:,j);
        for(int k=0; k<4; k++)
        {
            temp1[k] = m_birth[k][j];
        }
        matrizVector(F, 4, temp1, 4, out1, n);
        for(int k=0; k<4; k++)
        {
            m_birth[k][j] = Q[k][j] + out1[k];
        }
        calculateDataRange4(j, P_range, n);
        //P_birth(:,P_range) = G * Q * G' + F * P_birth(:,P_range) * F';
        for(int k=0; k<4; k++)
        {
            for(int i=0; i<4; i++) {
                temp[k][i] = P_birth[k][i];
            }
        }
        matrizxMatrizPB(G, 4, 4, Q, 4, 4, G, n, m);
        matrizxMatrizPB(G, 4, 4, G, 4, 4, G, n, m);
        matrizxMatrizPB(F, 4, 4, temp, 4, 4, temp, n, m);
        matrizxMatrizPB(temp, 4, 4, F, 4, 4, temp, n, m);
        for(int k=0; k<4; k++)
        {
            P_birth[k][j] = G[k][j] + temp[k][j];
        }
    }

    //Perform prediction for spawned targets using spawned velocities.
    for(int j = 0; j <numSpawnedTargets; j++)
    {

        calculate_Jacobian_F(dt, F, n);
        calculate_Jacobian_G(dt, G, n);
        //m_birth(:,j) = F * m_birth(:,j);
        for (int k = 0; k < 4; k++) {
            temp1[k] = m_spawn[k][j];
        }
        matrizVector(F, 4, temp1, 4, out1, n);
        for (int k = 0; k < 4; k++) {
            m_spawn[k][j] = Q[k][j] + out1[k];
        }
        calculateDataRange4(j, P_range, n);
        //P_birth(:,P_range) = G * Q * G' + F * P_birth(:,P_range) * F';
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 4; i++) {
                temp[k][i] = P_spawn[k][i];
            }
        }
        matrizxMatrizPB(G, 4, 4, Q, 4, 4, G, n, m);
        matrizxMatrizPB(G, 4, 4, G, 4, 4, G, n, m);
        matrizxMatrizPB(F, 4, 4, temp, 4, 4, temp, n, m);
        matrizxMatrizPB(temp, 4, 4, F, 4, 4, temp, n, m);
        for (int k = 0; k < 4; k++) {
            P_spawn[k][j] = G[k][j] + temp[k][j];
        }
    }

    if(VERBOSE == 1){
        for(int j = 0; j<numBirthedTargets; j++){
            //thisM = m_birth(:,j);
            for(int k=0; k<4; k++)
            {
                thisM[k] = m_birth[k][j];
            }
            cout << "Birthed target: " << j << thisM[0] << thisM[1] << thisM[2] << thisM[3] << "Weight "  << w_birth[j];
        }
        for(int j = 0; j<numSpawnedTargets; j++){
            //thisM = m_spawn(:,j);
            for(int k=0; k<4; k++)
            {
                thisM[k] = m_birth[k][j];
            }
            cout <<"Spawned target: " << j << thisM[0] << thisM[1] << thisM[2] << thisM[3] << "Weight "  << w_spawn[j];
        }
    }

}

void matrizxMatrizPB(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc){
    double sum=0;
    int i, j, k;

    if(ac != bf)
    {
        std::cout << "Different dimensions\n";
        exit(0);
    }
    nf = af;
    nc = bc;

    for(i=0;i<af;i++)
    {
        for(j=0;j<bf;j++)
        {
            for(k=0;k<bc;k++)
            {
                sum+=m1[i][k]*m2[k][j];
            }
            out[i][j]=sum;
            sum=0;
        }
    }

}
