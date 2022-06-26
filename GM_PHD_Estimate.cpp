//
// Created by david
//

/*
%GM_PHD_Estimate
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

%This file estimates the positions of targets tracked by the PHD filter.
%We need to extract the likely target positions from the PHD (i.e. we need to find the peaks of the PHD).
%This is actually fairly tricky. The naive approach is to pull out targets with the
%highest weights, but this is FAR from the best approach. A large covariance will pull down
%the peak size, and when targets are close together or have high covariances, there can be
%superposition effects which shift the peak.
*/

#include <iostream>
#include <cstring>
#include "GM_PHD_Estimate.h"
#include "GM_PHD_Initialisation.h"

using namespace std;

void GM_PHD_Estimate(void)
{
    //Previous definde variables
    extern double w_bar_k[9], weightThresholdToBeExtracted, m_bar_k[4][9], thisI, P_bar_k[4][36], X_k_history[NUM_F][NUM];
    extern int OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS, P_range[4], VERBOSE;

    //New global variables
    extern double X_k[4][3], X_k_P[4][12], X_k_w[3], i[3], thisP[4][4];
    extern int cX_k_P, cX_k, cX_k_w;

    int inum = 3, n;
    cX_k_P=0;
    cX_k=0;
    cX_k_w = 0;

    //This just implements the method in Vo&Ma, which is pulling out every target with a weight over
    //weightThresholdToBeExtracted (defined in GM_PHD_Initialisation). There is
    //the option of repeatedly printing out targets with rounded weights greater
    //than 1 (i.e. if two or more strong targets are mergde and the weight
    //rounds to 2/3/4/etc, display the target at that point multiple times when
    //VERBOSE is set to 1). This will NOT change filter performance as the
    //extracted state estimate is not fed back into the filter.
    cout<<"Step 6: Estimate target states" << endl;
    memset(X_k,0, 4*3*sizeof(double));
    memset(X_k_P,0, 4*12*sizeof(double));
    memset(X_k_w,0, sizeof(X_k_w));

    //OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS is set in GM_PHD_Initialisation
    if(OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS == 0)
    {
        for(int j=0; j<9; j++){
            if(w_bar_k[j] > weightThresholdToBeExtracted){
                i[j] = w_bar_k[j];
                inum++;
            }
            j++;
        }
        for(int k = 0; k<inum; k++)
        {
            for(int j = 0; j<inum; j++)
            {
                X_k[k][j] = m_bar_k[k][j];
                X_k_w[k] = w_bar_k[k];
            }
            cX_k++;
            cX_k_w++;

        }
        for(int j = 0; j<inum; j++)
        {
            thisI = i[j];
            calculateDataRange4(thisI,P_range, n);

            for(int k = 0; k<4; k++)
            {
                for (int i = 0; i < 4; i++)
                {
                    thisP[k][i] = P_bar_k[k][i];
                    X_k_P[k][cX_k_P + i] = thisP[k][i];
                }
            }
        }
    }
    else {
        //If a target has a rounded weight greater than 1, output it multiple
        //times. VERBOSE must be set to 1 to see the effects of this.
        for(int i = 0; i<9; i++) {
            for(int j = 0; j<1; j++) {
                X_k[i][cX_k + j] = m_bar_k[i][j];
                X_k_w[cX_k_w + i] = w_bar_k[i];
                calculateDataRange4(i, P_range, n);
                thisP[i][j] = P_bar_k[i][j];
                X_k_P[i][cX_k_P + j] = thisP[i][j];
            }
        }
    }

    if(VERBOSE == 1)
    {
        cout << "\t%d targets beleived to be valid: 3" << endl;
        for(int i = 0; i< 3; i ++)
        {
            calculateDataRange4(i, P_range, n);
            cout << "\t\tTarget " << i << " at " << X_k[0][i] << X_k[1][i] <<" P " << X_k_P[0][P_range[0]] << X_k_P[1][P_range[1]] << " W "
                         << X_k_w[i]<< endl;
        }
    }

    //Store history for plotting.
    //X_k_history = [X_k_history, X_k];
}
