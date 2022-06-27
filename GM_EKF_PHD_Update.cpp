//
// Created by david
//

#include <iostream>
#include "GM_EKF_PHD_Update.h"
#include "GM_PHD_Initialisation.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"

using namespace std;

void zerosU(int n, int m, double out[][NUM]);

void GM_EKF_PHD_Update(void) {

    //Previous define variables
    extern double w_k[NUM], m_k[NUM_F][NUM], P_k[NUM_F][NUM], Z[2][53], Pk_k_minus_1[NUM_F][NUM], prob_detection,
            wk_k_minus_1[NUM], mk_k_minus_1[NUM_F][NUM], m_j[4], mk_k_minus_1_before_prediction[NUM_F][NUM], dt, eta[NUM],
            S[4][NUM], K[4][NUM], P_k_k[4][NUM], simTarget1State[4], newP[4][4], thisZ[4], prevX[2], thisV[2], delta[4],
            m_new[4], P_new[4][4], thisEta[4], weight_tally, measZ[2], thisPos[4], x_sensor[3];;
    extern int numTargets_Jk_k_minus_1, P_range[4], VERBOSE, weightDataRange[2], numTargets_Jk, L, thisIndex, old_P_range[4],
            new_P_range[4], S_range[4], w_new, K_range[4], old_weight, new_weight, thisW;

    //New variables
    extern double measZ_cartesian[2];

    double out[1][NUM];
    int n;


    //This file performs a PHD filter update on the targets
    //This is basically a brute-force Kalman update of every target with every
    //measurement and creating a new target from the update results.

    cout << "Step 4: Performing update." << endl;

    //Set up matrices for post-update filter state
    zerosU(1, numTargets_Jk_k_minus_1 * 53 + numTargets_Jk_k_minus_1, out);
    for (int i = 0; i < numTargets_Jk_k_minus_1 * 53 + numTargets_Jk_k_minus_1; i++) {
        w_k[i] = out[0][i];
    }
    zerosU(4, numTargets_Jk_k_minus_1 * 53 + numTargets_Jk_k_minus_1, m_k);
    zerosU(4, 4 * (numTargets_Jk_k_minus_1 * 53 + numTargets_Jk_k_minus_1), P_k);

    //First we assume that we failed to detect all targets.
    //We scale all weights by probability of missed detection
    //We already did the prediction step for these so their position &
    //covariance will have been updated. What remains is to rescale their
    //weight.
    for (int j = 0; j < numTargets_Jk_k_minus_1; j++) {
        w_k[j] = (1 - prob_detection) * wk_k_minus_1[j];
        //m_k(:,j) = mk_k_minus_1(:,j);
        for (int k = 0; k < 4; k++) {
            m_k[k][j] = mk_k_minus_1[k][j];
        }
        calculateDataRange4(j, P_range, n);
        for (int k = 0; k < 4; k++) {
            newP[k][j] = Pk_k_minus_1[k][j];
        }
        for (int k = 0; k < 4; k++) {
            P_k[k][j] = newP[k][j];
        }
    }

    //Now we update all combinations of matching every observation with every target in the
    //map.
    //First we expand the observation to include velocity (by simple v = delta_x/delta_t for
    //delta_x = measured_position - position_before_prediction.
    //That is, how fast would the target have needed
    //to move to get from where it was to where it was seen now?
    L = 0;
    for(int zi = 0; zi<53; zi++)
    {
        L = L +
            1; //L is used to calculate an offset from previous updates. It maxes out at L = number_of_measurements. A more elegant solution would be to set L = zi but I retain this method for compatibility with Vo&Ma
        if (VERBOSE == 1){
            cout << "Matching targets against measurement " << zi << endl;
        }
        for(int j = 0; j<numTargets_Jk_k_minus_1; j++)
        {
            for (int k = 0; k < 4; k++) {
                m_j[k] = mk_k_minus_1[k][j];
            }
            //Augment the measurement of position with calculated velocity
            //This consists of only the observed position. But we need to extract the equivalent velocity observation
            for (int k = 0; k < 4; k++) {
                thisZ[k] = Z[k][j];
            }
            //Get the old (pre-prediction) position of the target
            for (int k = 0; k < 2; k++) {
                prevX[k] = mk_k_minus_1_before_prediction[k][j];
            }
            //velocity = dx / dt. Since Z and x are 2d, V = [Vx Vy]
            for (int k = 0; k < 2; k++) {
                thisV[k] = (thisZ[k]-prevX[k]) / dt;
            }
            //So we pretend to observe velocity as well as position
            for (int k = 2; k < 4; k++) {
                thisZ[k] = thisV[k];
            }

            //If this is the first time a target has been observed, it will have
            //no velocity stored.
            //Therefore it will be missing a prediction update, which will
            //impact both its position and its covariance.
            thisIndex = L * numTargets_Jk_k_minus_1 + j;

            calculateDataRange4(j, old_P_range, n);//Returns 4 columns
            for (int k = 0; k < 4; k++)
            {
                new_P_range[k] = 4 * L * numTargets_Jk_k_minus_1 + old_P_range[k];
            }
            calculateDataRange4(j, S_range, n);

            //Recalculate weight.
            //weightDataRange is used to control which dimensions we want to
            //reweight over. In my experience, reweighting over all four
            //produces unacceptably low weights most of the time, so we reweight
            //over 2D.
            /*
            w_new = prob_detection * wk_k_minus_1[j] * mvnpdf(thisZ[weightDataRange[0]], eta[weightDataRange[j]], S[weightDataRange[S_range[weightDataRange[0]]]]); //Hoping normpdf is the function I need
            w_k[thisIndex] = w_new;
             */

            //Update mean
            for(int k=0; k<4; k++)
            {
                delta[k] = thisZ[k] - eta[k];
            }
            calculateDataRange4(j, K_range, n);
            for(int k=0; k<4; k++)
            {
                m_new[k] = m_j[k] + K[k][j] * delta[k];

            }
            for(int k=0; k<thisIndex; k++)
            {
                m_k[k][j] = m_new[k];
            }

            //Update covariance
            for(int k=0; k<4; k++)
            {
                P_new[k][j] = P_k_k[k][j];
                P_k[k][j] = P_new[k][j];
            }

            if(VERBOSE == 1)
            {
                cout << "Observation: " << thisZ[0] << thisZ[1] << thisZ[2] << thisZ[3] << endl;
                for(int k=0; k<4; k++)
                {
                    thisEta[k] = eta[k];
                }
                cout << "Expected Obs: " << thisZ[0] << thisZ[1] << thisZ[2] << thisZ[3] << endl;

                cout << "Before Update: " << m_k[0] << m_k[1] << m_k[2] << m_k[3] << endl;
                cout << "After Update: " << m_k[0] << m_k[1] << m_k[2] << m_k[3] << endl;
                cout << "True State: " << simTarget1State[0] << simTarget1State[1] << simTarget1State[0] << simTarget1State[1]  << endl;
            }

            //Sum up weights for use in reweighting
            weight_tally = 0;
            for(int i = 0; i<numTargets_Jk_k_minus_1; i++)
            {
                thisIndex = L * numTargets_Jk_k_minus_1 + i;
                weight_tally = weight_tally + w_k[thisIndex];
            }

            //Recalculate weights
            if(VERBOSE == 1)
            {
                cout << "Calculating new weights for observation" << zi << endl;
            }

            for(int j = 0; j<numTargets_Jk_k_minus_1; j++)
            {
                old_weight = w_k[L * numTargets_Jk_k_minus_1 + j];
                inv_h(Z[0][zi], Z[1][zi], x_sensor[1], x_sensor[2], x_sensor[3], measZ_cartesian, n);
                new_weight = old_weight / (clutter_intensity(measZ_cartesian) + weight_tally); //Normalise

                w_k[L * numTargets_Jk_k_minus_1 + j] = new_weight;

                calculateDataRange4(j, S_range, n); //Data range for S, which is 4x4 matrix

            }


        }
        numTargets_Jk = L * numTargets_Jk_k_minus_1 + numTargets_Jk_k_minus_1;

        if(VERBOSE == 1)
        {
            for(int j = 0; j<numTargets_Jk; j++)
            {
                for(int k = 0; k<4; k++)
                {
                    thisPos[k] = m_k[k][j];
                }
                thisW = w_k[j];
                cout << "Target :" << j << thisPos[0] << thisPos[1] << thisPos[2] << thisPos[3] << "Weight" <<thisW;
            }
        }

    }
}

void zerosU(int n, int m, double out[][NUM])
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {
            out[i][j] = 0;
        }
    }
}