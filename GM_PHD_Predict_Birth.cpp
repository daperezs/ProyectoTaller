//
// Created by david on 06/06/2022.
//

/*
%GM_PHD_Predict_Birth
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
#include "GM_PHD_Predict_Birth.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"

using namespace std;

void GM_PHD_Predict_Birth(void){

    //Previos define variables
    extern double m_birth_before_prediction[NUM_F][NUM], m_birth[NUM_F][NUM], m_spawn[NUM_F][NUM], F[4][4],
            Q[4][4], P_birth[NUM_F][NUM], P_spawn[NUM_F][NUM], w_birth[NUM], w_spawn[NUM];
    extern int numBirthedTargets, numSpawnedTargets, i, VERBOSE;

    //New variables
    extern int P_range[4];
    extern double thisM[4];

    int n;
    double temp[4], out[4];

    //The means of these will be the position in Cartesian space where they are detected
    //The covariances & weights are calculated according to Vo&Ma
    cout<< "Step 1: Prediction for birthed and spawned targets." << endl;

    for(int i = 0; i < numBirthedTargets; i++)
    {
        m_birth_before_prediction[0][i] = m_birth[0][i];
        m_birth_before_prediction[1][i] = m_birth[1][i];
        m_birth_before_prediction[2][i] = m_birth[2][i];
        m_birth_before_prediction[3][i] = m_birth[3][i];
    }

    for(int i = 0; i < numSpawnedTargets; i++)
    {
        m_birth_before_prediction[0][numBirthedTargets+i] = m_spawn[0][i];
        m_birth_before_prediction[1][numBirthedTargets+i] = m_spawn[1][i];
        m_birth_before_prediction[2][numBirthedTargets+i] = m_spawn[2][i];
        m_birth_before_prediction[3][numBirthedTargets+i] = m_spawn[3][i];
    }


    //Perform prediction for birthed targets using birthed velocities.
    for(int j = 1; j <numBirthedTargets; j++)
    {
        int i = i + 1;
        //w_birth was already instantiated in GM_PHD_Create_Birth

        //m_birth(:,j) = F * m_birth(:,j);
        for(int k=0; k<4; k++)
        {
            temp[k] = m_birth[k][j];
        }
        matrizVector(F, 4, temp, 4, out, n);
        for(int k=0; k<4; k++)
        {
            m_birth[k][j] = Q[k][j] + out[k];
        }

        calculateDataRange4(j, P_range, n);

        //P_birth(:,P_range) = Q + F * P_birth(:,P_range) * F';
        for(int k=0; k<4; k++)
        {
            temp[k] = P_birth[k][P_range[k]];
        }
        matrizVector(F, 4, temp, 4, out, n);
        matrizVector(F, 4, out, 4, out, n);
        for(int k=0; k<4; k++)
        {
            P_birth[k][j] = Q[k][j] + out[k];
        }

    }

    //Perform prediction for spawned targets using spawned velocities.
    for(int j = 1; j <numSpawnedTargets; j++)
    {
        i = i + 1;
        //w_spawn was already instantiated in GM_PHD_Create_Birth
        //m_spawn(:,j) = F * m_spawn(:,j);
        for(int k=0; k<4; k++)
        {
            temp[k] = m_spawn[k][j];
        }
        matrizVector(F, 4, temp, 4, out, n);
        for(int k=0; k<4; k++)
        {
            m_spawn[k][j] = out[k];
        }

        calculateDataRange4(j, P_range, n);

        //P_spawn(:,P_range) = Q + F * P_spawn(:,P_range) *F';
        for(int k=0; k<4; k++)
        {
            temp[k] = P_spawn[k][P_range[k]];
        }
        matrizVector(F, 4, temp, 4, out, n);
        matrizVector(F, 4, out, 4, out, n);
        for(int k=0; k<4; k++)
        {
            P_spawn[k][j] = out[k];
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
