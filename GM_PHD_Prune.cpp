//
// Created by david
//

//GM_PHD_Prune
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

#include <iostream>
#include <cstring>
#include "GM_PHD_Prune.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"

using namespace std;

void max(double highWeights[], int n, double &maxW, int &j);
void sum(double w_k[], int iterateI, double w_bar_k_l);

void GM_PHD_Prune(void)
{

    // Previous define
    extern double w_k[NUM], T, P_k[NUM_F][NUM];
    extern int VERBOSE, P_range[4], mergeThresholdU, old_P_range[4];

    // New global variables
    extern double I[NUM], l, w_bar_k[9], m_bar_k[4][9], P_bar_k[4][36], highWeights[NUM], maxW, L[NUM], thisI,
            delta_m[4], mahal_dist, w_bar_k_l, m_bar_k_l[4], m_k[NUM_F][NUM], P_val[4][4], tmpP[4][4], P_bar_k_l[4][4],
            oldP[4][4];
    int j, n;

    int num = 0;

    //This file performs merging and pruning after the PHD update.
    //The PHD update creates a combinatorial explosion in the number of targets.
    //We prune the low-weight ones and merge the close-together ones. The weight
    //threshold T and distance threshold L that control this process are set in
    //GM_PHD_Initialisation.
    cout<< "Step 5: Prune and merge targets." << endl;

    ////Prune out the low-weight targets
    //Find targets with high enough weights
    int i = 0;
    while(i<3294){
        if(w_k[i] >= T){
            I[num] = w_k[i];
            num++;
        }
        i++;
    }

    if(VERBOSE == 1)
    {
        cout<< "The only tracks with high enough weights are:" << endl;
        cout<< I << endl;
    }

    //// Merge the close-together targets
    l = 0; //Counts number of features
    memset(w_bar_k,0, sizeof(w_bar_k));
    memset(m_bar_k,0, 4*9*sizeof(double));
    memset(P_bar_k,0, 4*36*sizeof(double));
    while(!isempty(I, num)){ //We delete from I as we merge
        l = l + 1;
        //Find j, which is i corresponding to highest w_k for all i in I
        for(int k = 0; k<num; k++)
            highWeights[k] = w_k[k];

        max(highWeights, num, maxW, j);
        //In case of two targets with equal weight
        //j is an index of highWeights (i.e. a position in I)
        //We want the index in w_k
        j = I[j];

        //Find all points with Mahalanobis distance less than U from point
        //m_k(j)
        memset(L,0, sizeof(L));//A vector of indexes of w_k
        for(int iterateI = 0; iterateI<num; iterateI++)
        {
            thisI = I[iterateI];

            for(int k = 0; k<4; k++){
                //delta_m[k] = m_k(:,thisI) - m_k(:,j);
            }
            calculateDataRange4(thisI, P_range, n);
            double aux[4];
            for(int k = 0; k<4; k++){
                for(int j = 0; j<4; j++){
                    aux[k] = aux[k] + P_k[k][j];
                }
            }
            for(int k = 0; k<4; k++)
                mahal_dist = delta_m[k] * (aux[k] / delta_m[k]); //Invert covariance via left division

            if (mahal_dist <= mergeThresholdU){
                L[iterateI] = thisI;
            }

            if(VERBOSE == 1){
                cout<< "\tMerging target %d with these targets:" << j << endl;
                cout<< L << endl;
            }

            //The new weight is the sum of the old weights
            sum(w_k, iterateI, w_bar_k_l);

            //The new mean is the weighted average of the merged means
            memset(m_bar_k_l,0, sizeof(m_bar_k_l));
            for(int i = 0; i < iterateI; i++){
                thisI = L[i];
                //m_bar_k_l[i] = m_bar_k_l + 1 / w_bar_k_l *  (w_k[thisI] * m_k[thisI]);
            }

            //Calculating covariance P_bar_k is a bit trickier
            memset(P_val,0, 4*4*sizeof(double));
            for(int i = 0; i< iterateI; i++){
                thisI = L[i];
                //delta_m = m_bar_k_l - m_k(:,thisI);
                calculateDataRange4(thisI, P_range, n);
                for(int k = 0; k<4; k++) {
                    for (int j = 0; j < 4; j++) {
                        P_val[k][j] = P_val[k][j] + w_k[k] * (P_k[k][j] + delta_m[k] * delta_m[k]);
                        tmpP[k][j] = P_k[k][j];
                        P_bar_k_l[k][j] = P_val[k][j] / w_bar_k_l;
                    }
                }
            }

            calculateDataRange4(j, old_P_range, n);
            for(int k = 0; k<4; k++) {
                for (int j = 0; j < 4; j++) {
                    oldP[k][j] = P_k[k][j];
                }
            }

            //Now delete the elements in L from I
            for(int i = 0; i< iterateI; i++)
            {
                //iToRemove = find(I == L(i));
                //I[iToRemove] = [];
            }


            //Append the new values to the lists
            for(int k = 0; k<4; k++)
            {
                w_bar_k[k] = w_bar_k_l;

                m_bar_k[0][k] = m_bar_k_l[k];
                m_bar_k[1][k] = m_bar_k_l[k];
                m_bar_k[2][k] = m_bar_k_l[k];
                m_bar_k[3][k] = m_bar_k_l[k];

                P_bar_k[0][k] = P_bar_k_l[0][k];
                P_bar_k[1][k] = P_bar_k_l[1][k];
                P_bar_k[2][k] = P_bar_k_l[2][k];
                P_bar_k[3][k] = P_bar_k_l[3][k];
            }

        }

    }


}

void max(double highWeights[], int n, double &maxW, int &j){
    for(int k = 0; k<n; k++)
    {
        if(highWeights[k] > maxW)
        {
            maxW = highWeights[k];
            j = k+1;
        }
    }
}

void sum(double w_k[], int iterateI, double w_bar_k_l)
{
    for(int i = 0; i < iterateI; i++)
    {
        w_bar_k_l = w_bar_k_l + w_k[iterateI];
    }
}
