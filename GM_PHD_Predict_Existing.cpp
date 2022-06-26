//
// Created by david on 17/06/2022.
//

//GM_PHD_Predict_Existing
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

#include <iostream>
#include <cstring>
#include "GM_PHD_Predict_Existing.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"

using namespace std;


void GM_PHD_Predict_Existing(void)
{
    //Previos define variables
    extern double mk_k_minus_1_before_prediction[NUM_F][NUM], mk_minus_1[NUM_F][NUM], wk_minus_1[NUM], prob_survival,
            F[4][4], Q[4][4], Pk_minus_1[NUM_F][NUM], wk_k_minus_1[NUM], mk_k_minus_1[NUM_F][NUM], Pk_k_minus_1[NUM_F][NUM],
            w_birth[NUM], m_birth[NUM_F][NUM], P_birth[NUM_F][NUM], w_spawn[NUM], m_spawn[NUM_F][NUM], P_spawn[NUM_F][NUM],
            m_birth_before_prediction[NUM_F][NUM];
    extern int P_range[4], VERBOSE, numTargets_Jk_k_minus_1, numTargets_Jk_minus_1,  numBirthedTargets, numSpawnedTargets;

    //New variables
    extern double P_i[4][4], prevState[4], newState[4];
    int n;
    double temp[4], out[4];

    //This file performs prediction for existing targets
    cout<< "Step 2: Prediction for existing targets." << endl;

    memcpy(mk_k_minus_1_before_prediction, mk_minus_1, 4 * sizeof(double));

    for(int j = 0; j<9; j++)
    {
        wk_minus_1[j] = prob_survival * wk_minus_1[j];
        //mk_minus_1(:,j) = F * mk_minus_1(:,j); //Assume constant velocity.
        for(int k=0; k<4; k++)
        {
            temp[k] = mk_minus_1[k][j];
        }
        matrizVector(F, 4, temp, 4, out, n);
        for(int k=0; k<4; k++)
        {
            mk_minus_1[k][j] = out[k];
        }

        calculateDataRange4(j, P_range, n);
        //P_i = Q + F * Pk_minus_1(:,P_range) * F;
        for(int k=0; k<4; k++)
        {
            temp[k] = Pk_minus_1[k][P_range[k]];
        }
        matrizVector(F, 4, temp, 4, out, n);
        matrizVector(F, 4, out, 4, out, n);
        for(int k=0; k<4; k++)
        {
            P_i[k][j] = Q[k][j] + out[k];
        }

        for(int k=0; k<4; k++)
        {
            prevState[k] = mk_k_minus_1_before_prediction[k][j];
            newState[k] = mk_minus_1[k][j];
        }

        for(int k=0; k<4; k++)
        {
            Pk_minus_1[k][P_range[j]] =  P_i[k][j];
        }

        if(VERBOSE == 1)
        {
            cout<<"\t\tExisting target %d. Previously at %3.4f %3.4f, now at %3.4f %3.4f." <<j << prevState[1]<< prevState[2]<< newState[1]<< newState[2]<< endl;

            cout<< "\t\tP was %3.4f %3.4f, NOW %3.4f %3.4f" << Pk_minus_1[1][P_range[1]]<< Pk_minus_1[2][P_range[2]]<< P_i[1][1]<< P_i[2][2]<< endl;

        }
    }

    //Now we combine the birthed targets with the existing ones.
    //Append newly birthed targets (in m_k_minus_1) to back of old ones

    int size1 = 0;
    for(int i= 0; i< 9; i++){
        wk_k_minus_1[i]= wk_minus_1[i];
        size1 = i+1;
    }
    int size11 = 0;
    for(int i= 0; i< 9; i++){
        wk_k_minus_1[size1+i]= w_birth[i];
        size11 = size1 + i +1;
    }
    for(int i= 0; i< 9; i++){
        wk_k_minus_1[size11+i]= w_spawn[i];
    }

    int size2 = 0;
    for(int i= 0; i<9; i++){
        mk_k_minus_1[0][i] = mk_minus_1[0][i];
        mk_k_minus_1[1][i] = mk_minus_1[1][i];
        mk_k_minus_1[2][i] = mk_minus_1[2][i];
        mk_k_minus_1[3][i] = mk_minus_1[3][i];
        size2 = i +1;
    }
    int size3 = 0;
    for(int i= 0; i< 9; i++){
        mk_k_minus_1[0][size2+i] = m_birth[0][i];
        mk_k_minus_1[1][size2+i] = m_birth[1][i];
        mk_k_minus_1[2][size2+i] = m_birth[2][i];
        mk_k_minus_1[3][size2+i] = m_birth[3][i];
        size3=size2+i+1;
    }
    for(int i= 0; i< 9; i++){
        mk_k_minus_1[0][size3+i] = m_spawn[0][i];
        mk_k_minus_1[1][size3+i] = m_spawn[1][i];
        mk_k_minus_1[2][size3+i] = m_spawn[2][i];
        mk_k_minus_1[3][size3+i] = m_spawn[3][i];
    }

    size2 = 0;
    for(int i= 0; i<9; i++){
        Pk_k_minus_1[0][i] = Pk_minus_1[0][i];
        Pk_k_minus_1[1][i] = Pk_minus_1[1][i];
        Pk_k_minus_1[2][i] = Pk_minus_1[2][i];
        Pk_k_minus_1[3][i] = Pk_minus_1[3][i];
        size2 = i +1;
    }
    size3 = 0;
    for(int i= 0; i< 9; i++){
        Pk_k_minus_1[0][size2+i] = P_birth[0][i];
        Pk_k_minus_1[1][size2+i] = P_birth[1][i];
        Pk_k_minus_1[2][size2+i] = P_birth[2][i];
        Pk_k_minus_1[3][size2+i] = P_birth[3][i];
        size3=size2+i+1;
    }
    for(int i= 0; i< 9; i++){
        Pk_k_minus_1[0][size3+i] = P_spawn[0][i];
        Pk_k_minus_1[1][size3+i] = P_spawn[1][i];
        Pk_k_minus_1[2][size3+i] = P_spawn[2][i];
        Pk_k_minus_1[3][size3+i] = P_spawn[3][i];
    }

    numTargets_Jk_k_minus_1 = numTargets_Jk_minus_1 + numBirthedTargets + numSpawnedTargets;
    //Create a backup to allow for augmenting the measurement in the update
    for(int i = 0; i<9; i++){
        mk_k_minus_1_before_prediction[0][i] = m_birth_before_prediction[0][i];
        mk_k_minus_1_before_prediction[1][i] = m_birth_before_prediction[1][i];
        mk_k_minus_1_before_prediction[2][i] = m_birth_before_prediction[2][i];
        mk_k_minus_1_before_prediction[3][i] = m_birth_before_prediction[3][i];
    }
    //mk_k_minus_1_before_prediction = [mk_k_minus_1_before_prediction, m_birth_before_prediction]; //m_birth_before_prediction also contains the spawned targets before prediction

    if(VERBOSE == 1)
    {
        cout << "\tPerformed prediction for %d birthed and existing targets in total." << numTargets_Jk_k_minus_1 << endl;
    }
}
