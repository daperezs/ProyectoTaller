//
// Created by david on 04/05/2022.
//

#ifndef PROYECTO_GM_VARGLOBALES_H
#define PROYECTO_GM_VARGLOBALES_H

#define NUM_F 10
#define NUM 30000

//GM_PHD_Initialisation
int VERBOSE, KNOWN_TARGET, PLOT_ALL_MEASUREMENTS, OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS,
        CALCULATE_OSPA_METRIC, USE_EKF, DRAW_VELOCITIES, addVelocityForNewTargets, addStaticNewTargets,
        weightDataRange[2], numBirthedTargets, numSpawnedTargets, numTargets_Jk_k_minus_1, numTargets_Jk,
        numTargets_J_pruned, numTargets_Jk_minus_1, mergeThresholdU, maxGaussiansJ, xrange[2], yrange[2],
        z_cartesian[2],order_p, cutoff_c, c1Birth, c2Birth, c3Birth;
double temp[NUM], MAX_V, m_birth_before_prediction[NUM_F][NUM], prob_survival,
        mk_k_minus_1_before_prediction[NUM_F][NUM], wk_k_minus_1[NUM], mk_k_minus_1[NUM_F][NUM],
        Pk_k_minus_1[NUM_F][NUM], w_k[NUM], m_k[NUM_F][NUM], P_k[NUM_F][NUM], T, weightThresholdToBeExtracted,
        wk_minus_1[NUM], mk_minus_1[NUM_F][NUM], Pk_minus_1[NUM_F][NUM], X_k_history[NUM_F][NUM], w_birth[NUM],
        m_birth[NUM_F][NUM], P_birth[NUM_F][NUM], w_spawn[NUM], m_spawn[NUM_F][NUM], P_spawn[NUM_F][NUM],
        V,lambda_c, metric_history[NUM],dt, F[4][4], sigma_v, Q[4][4],
        birth_mean1[4], birth_mean2[4], covariance_birth[4][4], covariance_spawn[4][4], prob_detection, H[2][4],
        H2[4][4], R2[4][4], sigma_r, k;

//GM_PHD_Simulate_Initialise
int nClutter, endTime,c1History, c2History, c3History;
double noiseScaler, simTarget1Start[4], simTarget2Start[4], simTarget3Start[4],simTarget1End[2],
        simTarget2End[2], simTarget3End[2], simTarget1Vel[2], simTarget2Vel[2], simTarget3Vel[2],
        simTarget1History[NUM_F][NUM], simTarget2History[NUM_F][NUM], simTarget3History[NUM_F][NUM],
        simTarget1State[4], simTarget2State[4], simTarget3State[4], simTarget3SpawnTime;

//GM_EKF_PHD_Simulate_Initialise
double x_sensor[3];

//GM_EKF_PHD_Initialise_Jacobians
double sigma_theta, upsilon_x, upsilon_y, upsilon_vx, upsilon_vy, U[4][4];

//GM_PHD_Simulate_Measurements
double clutterX, clutterY, detect1, detect2, detect3, measX1, measY1, measX2, measY2, measX3,
        measY3, randn, Z[2][53], zTrue[2][3];

//GM_EKD_PHD_Simulate_Measurements
double clutterRTheta[2], measR1, measTheta1, simMeas1[2], measR2, measTheta2, simMeas2[2], measR3, measTheta3,
        simMeas3[2], clutter[NUM_F][NUM];

//GM_PHD_Predict_Birth
int P_range[4];
double thisM[4];

//GM_PHD_Predict_Existing
double P_i[4][4], prevState[4], newState[4];

//GM_PHD_Construct_Update_Components
double eta[NUM], S[4][NUM], K[4][NUM], P_k_k[4][NUM], m_j[4], eta_j[4], PHt[4][4], S_j[4][4], SChol[4][4], SCholInv[4][4],
            K_j[4][4], W1[4][4], P_j[4][4];
int ceta, cS, cK, cP_k_k;


//GM_PHD_Update
double newP[4][4], thisZ[4], prevX[2], thisV[2], delta[4], m_new[4], P_new[4][4], thisEta[4], weight_tally,
            measZ[2], thisPos[4];
int L, thisIndex, old_P_range[4], new_P_range[4], S_range[4], w_new, K_range[4], old_weight, new_weight, thisW;

//GM_EKF_PHD_Predict_Birth
double G[4][4];

//GM_EKF_PHD_Construct_Update_Components
double xR, yR, hR, R_polar[2][2], r, theta, R_cartesian[2][2], R_velocity[2][2], R[4][4], HPHt[4][4];

//GM_EKF_PHD_Update
double measZ_cartesian[2];

//GM_PHD_Prune
double I[NUM], l, w_bar_k[9], m_bar_k[4][9], P_bar_k[4][36], highWeights[NUM], maxW, thisI,
            delta_m[4], mahal_dist, w_bar_k_l, m_bar_k_l[4], P_val[4][4], tmpP[4][4], P_bar_k_l[4][4],
            oldP[4][4];

// GM_PHD_Estimate
double X_k[4][3], X_k_P[4][12], X_k_w[3], i[3], thisP[4][4];
int cX_k_P, cX_k, cX_k_w;

//GM_PHD_Create_Birth
double thisMeasRowRange, prevMeasRowRange;

//GM_PHD_Calculate_Performance_Metric
double X[4][3], Y[4][3], metric;
int cY;

#endif //PROYECTO_GM_VARGLOBALES_H
