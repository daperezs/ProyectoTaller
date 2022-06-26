//
// Created by david
//

//GM_PHD_Calculate_Performance_Metric
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

#include <cstring>
#include "GM_PHD_Calculate_Performance_Metric.h"

void GM_PHD_Calculate_Performance_Metric(void)
{
    //Previos define variables
    extern int CALCULATE_OSPA_METRIC, order_p, cutoff_c;
    extern double X_k[4][3], simTarget1History[NUM_F][NUM], simTarget2History[NUM_F][NUM], simTarget3History[NUM_F][NUM],
            k, simTarget3SpawnTime;

    //New global variables
    extern double X[4][3], Y[4][3], metric, metric_history[100];
    extern int cY;

    cY=0;

    //Schuhmacher, D.; Ba-Tuong Vo; Ba-Ngu Vo, "A Consistent Metric for Performance Evaluation of Multi-Object Filters," Signal Processing, IEEE Transactions on , vol.56, no.8, pp.3447,3457, Aug. 2008
    //cutoff_c and order_p are tuning parameters set in GM_PHD_Initialisation
    if(CALCULATE_OSPA_METRIC == 1)
    {
        memcpy(X,X_k, 4*3*sizeof(double));
        for(int i = 0; i<k; i++){
            for(int j = 0; j<k; j++){
                Y[i][j] = simTarget1History[i][j];
                cY++;
            }
        }
        for(int i = k; i<2*k; i++){
            for(int j = k; j<2*k; j++){
                Y[i][j] = simTarget2History[i][j];
                cY++;
            }
        }
    }
    if(k >= simTarget3SpawnTime)
    {
        for(int i = 0; i<k-simTarget3SpawnTime+1; i++){
            for(int j = 0; j<k-simTarget3SpawnTime+1; j++){
                Y[i][cY + j] = simTarget3History[i][j];
                cY++;
            }
        }
    }


    //metric = ospa_dist(X, Y, cutoff_c, order_p); //Use Ba-Ngu Vo's implementation
    //metric_history = [metric_history, metric];
    //This uses the Optimal Subpattern Assignment (OSPA) metric proposed by

}