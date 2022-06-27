//
// Created by david
//

#include <cstring>
#include <cmath>
#include "GM_PHD_Initialisation.h"
#include "unifpdf_2d.h"

//GM_PHD_Initialisation
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

//This file initialises most of the variables that we use for the filter.

//If you want to use this GM-PHD filter for your own problem, you will need
//to replace a lot of this script with your own code.

void GM_PHD_Initialisation(void)
{
    // declaramos las variables globales
    extern int VERBOSE, KNOWN_TARGET, PLOT_ALL_MEASUREMENTS, OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS,
            CALCULATE_OSPA_METRIC, USE_EKF, DRAW_VELOCITIES, addVelocityForNewTargets, addStaticNewTargets,
            weightDataRange[2], numBirthedTargets, numSpawnedTargets, numTargets_Jk_k_minus_1, numTargets_Jk,
            numTargets_J_pruned, numTargets_Jk_minus_1, mergeThresholdU, maxGaussiansJ, xrange[2], yrange[2],
            z_cartesian[2],order_p, cutoff_c, c1Birth, c2Birth, c3Birth;
    extern double temp[NUM], MAX_V, m_birth_before_prediction[NUM_F][NUM], prob_survival,
            mk_k_minus_1_before_prediction[NUM_F][NUM], wk_k_minus_1[NUM], mk_k_minus_1[NUM_F][NUM],
            Pk_k_minus_1[NUM_F][NUM], w_k[NUM], m_k[NUM_F][NUM], P_k[NUM_F][NUM], T, weightThresholdToBeExtracted,
            wk_minus_1[NUM], mk_minus_1[NUM_F][NUM], Pk_minus_1[NUM_F][NUM], X_k_history[NUM_F][NUM], w_birth[NUM],
            m_birth[NUM_F][NUM], P_birth[NUM_F][NUM], w_spawn[NUM], m_spawn[NUM_F][NUM], P_spawn[NUM_F][NUM],
            V,lambda_c, metric_history[NUM],dt, F[4][4], sigma_v, Q[4][4], k,
            birth_mean1[4], birth_mean2[4], covariance_birth[4][4], covariance_spawn[4][4], prob_detection, H[2][4],
            H2[4][4], R2[4][4], sigma_r;

    int m;

    // These variables control some aspect of filter function and performance
    VERBOSE = 0; //Set to 1 for much more text output. This can be used to give more information about the state of the
    // filter, but slows execution and makes the code less neat by requiring disp() statements everywhere.
    KNOWN_TARGET = 1; //Set to 1 to have the first two targets already known to the filter on initialisation. Otherwise
    // tracks for them will need to be instantiated.
    PLOT_ALL_MEASUREMENTS = 0; //Set to 1 to maintain the plot of the full measurement history, including clutter and
    // error ellipses. Set to 0 to just plot the most recent measurement.
    OUTPUT_MULTIPLE_HIGH_WEIGHT_TARGETS = 0; //Set to 1 to implement Vo & Ma's state extraction where targets with
    // rounded weight greater than 1 are output multiple times. This does NOT change filter processing as extracted
    // targets are not fed back into the filter, it is merely for display. VERBOSE must be set to 1 to see the effects
    // of this. Personally I think this is a bit redundant but I've included it to match the Vo & Ma algorithm.
    CALCULATE_OSPA_METRIC = 1; //Set to 1 to calculate the OSPA performance metric for each step. This is not essential
    // for the GM-PHD filter but provides a rough indication of how well the filter is performing at any given timestep.
    USE_EKF = 0; //Set to 1 to use extended Kalman filter. Set to 0 to use linear KF.
    DRAW_VELOCITIES = 0; //Set to 1 to draw velocity arrows on the output plot.


    //Target initialisation: when we add a new target, we can use a two-step
    //initialisation where every target is added after two observations, so it
    //is added with a position (the second observation position) and a velocity
    //(the velocity from the first to the second observation). Or we can add
    //targets as static objects after one observation, and let the velocity be
    //filled in by subsequent observations. The big difference is the first
    //prediction step - with no velocity initially they will be expected to stay
    //in the same location, as opposed to move in a straight line.
    //We can add both types of target, but using a two-step initialisation: so a
    //after the second observation, a static target is added and a target with
    //velocity.
    addVelocityForNewTargets = 0;
    addStaticNewTargets = 1;

    //When recalculating the weights of targets, we consider the match between
    //observed and predicted position. We could use the observed and predicted
    //velocity as well, but I haven't gotten around to implementing that. So for
    //reweighting, we only use the first two values of the state.
    // weightDataRange = 1:2;
    weightDataRange[0] = 1;
    weightDataRange[1] = 1;

    memset(temp,0, sizeof(temp));

    //Filtering velocity: We throw away observed targets with velocities greater
    //than this.
    MAX_V = 100;

    //Data structures and variables
    //Step 0: Initialisation
    k = 0; //Time step

    //Step 1: Prediction for birthed targets
    numBirthedTargets = 0;
    numSpawnedTargets = 0;

    c1Birth=0;
    c2Birth=0;
    c3Birth=0;
    memset(m_birth_before_prediction,0, NUM_F*NUM*sizeof(double)); //We store the birth/spawn position
    // from before the prediction. This is used in augmenting the measurement vector to calculate velocity, for the
    // update.

    //Step 2: Prediction for existing targets
    //We store the existing position from before the prediction. This is used in augmenting measurement vector to
    // calculate velocity, for the update.
    memset(mk_k_minus_1_before_prediction,0, NUM_F*NUM*sizeof(double));//Used in augmenting measurement
    // vector to calculate velocity, for update.
    numTargets_Jk_k_minus_1 = 0;//Number of targets given previous. J_k|k-1. Set at end of Step 2 (prediction of
    // existing targets)
    prob_survival = 0.99; //Probability of target survival. Used in GM_PHD_Predict_Existing for weight calculation

    //These are used for storing the state of existing targets after prediction.
    memset(wk_k_minus_1,0, sizeof(wk_k_minus_1));//Weights of gaussians, previous, predicted. w_k|k-1.
    memset(mk_k_minus_1,0, NUM_F*NUM*sizeof(double));//Means of gaussians, previous, predicted. m_k|k-1
    memset(Pk_k_minus_1,0, NUM_F*NUM*sizeof(double));//Covariances of gaussians, previous, predicted. P_k|k-1

    //Step 3: Construct update component.
    //There is nothing we need to initialise here because the arrays are initialised
    //each iteration. (This is actually true for a lot of the arrays being initialised in this script arrays, but not
    //all, and rather than try to separate them it's easiest to initialise nearly all of them here).

    //Step 4: Update
             //The posterior weight, mean and covariance after the update
    memset(w_k,0, sizeof(w_k));//Weights of gaussians, updated. w_k|k
    memset(m_k,0, NUM_F*NUM*sizeof(double));//Means of gaussians, updated. m_k|k
    memset(P_k,0, NUM_F*NUM*sizeof(double));//covariances of gaussians, updated. P_k|k
    numTargets_Jk = 0; //Number of targets after update. J_k. Set at end of step 4 (update)


    //Step 5: Prune
    numTargets_J_pruned = 0;//Number of targets after pruning
    numTargets_Jk_minus_1 = 0; //Number of targets, previous. J_k-1. Set in end of GM_PHD_Prune
    //Merge and prune constants
    //These numbers have a HUGE impact on the performance and unfortunately need to be
    //manually tuned if we make changes.
    //These particular values come from Vo&Ma
    T = 10e-5;//Weight threshold. Value the weight needs to be above to be considered a target rather than be deleted
    // immediately.
    mergeThresholdU = 4; //Merge threshold. Points with Mahalanobis distance of less than this between them will be merged.
    weightThresholdToBeExtracted = 0.5;//Value the weight needs to be above to be considered a 'real' target.
    maxGaussiansJ = 100;//Maximum number of Gaussians after pruning. NOT USED in this implementation.

    //The previous iteration's mean/weight/covariance. Set in GM_PHD_Prune
    //after pruning. Used as input for the prediction step.
    memset(wk_minus_1,0, sizeof(double));//Weights from previous iteration
    memset(mk_minus_1,0, NUM_F*NUM*sizeof(double));//Means from previous iteration
    memset(Pk_minus_1,0, NUM_F*NUM*sizeof(double));//Covariances from previous iteration

    //Step 6: Estimate/extract states
    //We store the history of all points X_k (extracted states) for plotting
    //purposes. This is updated in the end of GM_PHD_Estimate
    memset(X_k_history,0, NUM_F*NUM*sizeof(double));

    //Step 7: Create birthed/spawned targets to append to list next iteration
    //These are set in GM_PHD_Create_Birth
    //Births occur at fixed positions; spawns occur at existing targets
    memset(w_birth,0, sizeof(double)); //New births' weights
    memset(m_birth,0, NUM_F*NUM*sizeof(double)); //New births' means
    memset(P_birth,0, NUM_F*NUM*sizeof(double)); //New births' covariances
    memset(w_spawn,0, sizeof(double)); //New spawns' weights
    memset(m_spawn,0, NUM_F*NUM*sizeof(double)); //New spawns' means
    memset(P_spawn,0, NUM_F*NUM*sizeof(double)); //New spawns' covariances


    //Step Sim: Generate simulated measurements
    //Detected clutter is a Poisson RFS
    xrange[0] = -1000; //X range of measurements
    xrange[1] = 1000; //X range of measurements
    yrange[0] = -1000; //Y range of measurements
    yrange[1] = 1000; //Y range of measurementsts
    V = 4 * 10e6; //Volume of surveillance region
    lambda_c = 12.5 * 10e-6; //average clutter returns per unit volume (50 clutter returns over the region)
    //clutter_intensity = lambda_c * V * unifpdf_2d(xrange, yrange, z_cartesian);//Generate clutter function. There are
    // caveats to its use for clutter outside of xrange or yrange - see the comments in unifpdf_2d.m

    //Step Metric: Calculate performance metric
    //The order_p and cutoff_c control filter performance; read the paper by
    //Schuhmacher et al to understand these in greater detail. The parameters
    //given here by default are not particularly well tuned for this problem but
    //they work alright.
    order_p = 1;//The order determines how punishing the metric calculation is to larger errors; as p increases,
    // outliers are more heavily penailsed
    cutoff_c = 200;//Cutoff determines the maximum error for a point.
    memset(metric_history,0, sizeof(metric_history)); //History of the OSPA performance metric

    //MODELS for prediction and observation
    //Prediction models - used in steps 1 & 2 for prediction
    //eye(2, I2);//2x2 identify matrix, used to construct matrices
    //memset(Z2,0, NUM*NUM*sizeof(Z2));//2x2 zero matrix, used to construct matrices
    dt = 1; //One-second sampling period
    //F = [ [I2, dt*I2]; [Z2 I2] ]; //State transition matrix (motion model)
    memset(F,0, 4*4*sizeof(double));
    for(int i=0; i<4; i++) F[i][i]=1;
    F[0][2]= F[1][3] = dt;
    sigma_v = 5; //Standard deviation of process noise is 5 m/(s^2)
    //Q = sigma_v^2 * [ [1/4*dt^4*I2, 1/2*dt^3*I2]; [1/2*dt^3* I2, dt^2*I2] ]; //Process noise covariance, given in Vo&Ma.
    memset(Q,0, 4*4*sizeof(double));
    Q[0][0] = Q[1][1] = 1.0/4.0*sigma_v*sigma_v*pow(dt,4);
    Q[2][2] = Q[3][3] = sigma_v*sigma_v*pow(dt,2);
    Q[0][2] = Q[1][3] = 1.0/2.0*sigma_v*sigma_v*pow(dt,3);
    Q[2][0] = Q[3][1] = 1.0/2.0*sigma_v*sigma_v*pow(dt,3);

    //Birth model. This is a Poisson random finite set.
    //We only use the first two elements as the initial velocities are unknown
    //and it is not specified by Vo&Ma whether they are used (I am fairly sure
    //they aren't, or they would have said otherwise)
    //Birth and spawn models
    //birth_mean1 = [250, 250, 0, 0];//Used in birth_intensity function
    //birth_mean2 = [-250, -250, 0, 0];//Used in birth_intensity function
    birth_mean1[0] = birth_mean1[1] = 250;
    birth_mean1[2] = birth_mean1[3] = 0;
    birth_mean2[0] = birth_mean2[1] = -250;
    birth_mean2[2] = birth_mean2[3] = 0;
    //covariance_birth = diag([100, 100, 25, 25]');%Used in birth_intensity function
    //covariance_spawn = diag([100, 100, 400,400]');%Used in spawn_intensity function
    memset(covariance_birth,0, 4*4*sizeof(double)); //Used in birth_intensity function
    covariance_birth[0][0] = covariance_birth[1][1] = 100;
    covariance_birth[2][2] = covariance_birth[3][3] = 25;
    memset(covariance_spawn,0, 4*4*sizeof(double)); //Used in spawn_intensity function
    covariance_spawn[0][0] = covariance_spawn[1][1] = 100;
    covariance_spawn[2][2] = covariance_spawn[3][3] = 400;
    //covariance_spawn = max(covariance_spawn, 10^-6);//Used in spawn_intensity function
    max(covariance_spawn,4,1e-6,covariance_spawn,m);
    //birth_intensity = @(x) (0.1 * mvnpdf(x(1:2), birth_mean1(1:2), covariance_birth(1:2,1:2)) + 0.1 * mvnpdf(x(1:2),birth_mean2(1:2), covariance_birth(1:2,1:2)));//Generate birth weight. This only takes into account the position, not the velocity, as Vo&Ma don't say if they use velocity and I assume that they don't. Taken from page 8 of their paper.
    //spawn_intensity = @(x, targetState) 0.05 * mvnpdf(x, targetState, covariance_spawn);//Spawn weight, from page 8 of Vo&Ma.

    prob_detection = 0.98; //Probability of target detection. Used in recalculating weights in GM_PHD_Update

    //Detection models for the linear Kalman filter. The extended Kalman filter
    //uses different models.
    if(USE_EKF == 0)
    {
        //H = [I2, Z2];//Observation matrix for position. Not used, but if you wanted to cut back to just tracking position, might be useful.
        memset(H,0, 2*4*sizeof(double));
        H[0][0] = H[1][1] = 1;
        //H2 = eye(4);//Observation matrix for position and velocity. This is the one we actually use, in GM_PHD_Construct_Update_Components
        memset(H2,0, 4*4*sizeof(double));
        H2[0][0] = H2[1][1] = H2[2][2] = H2[3][3] = 1;
        sigma_r = 10; //Standard deviation of measurement noise is 10m. Used in creating R matrix (below)
        //R = sigma_r ^ 2 * I2;//Sensor noise covariance. used in R2 (below)
        //R2 = [[R, Z2];
        //[Z2, R] ];//Measurement covariance, expanded to both position & velocity. Used in GM_PHD_Construct_Update_Components. NOTE: This assumes that speed measurements have the same covariance as position measurements. I have no mathematical justification for this.
        memset(R2,0, 4*4*sizeof(double));
        R2[0][0] = R2[1][1] = R2[2][2] = R2[3][3] = sigma_r*sigma_r;
    }
}

//Utility functions
//Used to calculate the indices of two-dimensional target j in a long list of two-dimensional targets
void calculateDataRange2(int j, int out[], int &n)
{
    n=2;
    out[0] = 2*j-1;
    out[1] = 2*j;
}

//Used to calculate the indices of four-dimensional target j in a long list of four-dimensional targets
void calculateDataRange4(int j, int out[], int &n)
{
    n=4;
    out[0] = 4*(j-1) + 1;
    out[1] = 4*(j-1) + 2;
    out[2] = 4*(j-1) + 3;
    out[3] = 4*(j-1) + 4;
}

//Generate clutter function. There are caveats to its use for clutter outside of xrange or yrange - see the comments in unifpdf_2d.m
double clutter_intensity(double z[])
{
    extern int xrange[2], yrange[2];
    extern double lambda_c, V;

    return(lambda_c * V * unifpdf_2d(xrange, yrange, z));
}

void max(double a[][4], int n, double num, double out[][4], int &m)
{
    m=n;
    memcpy(out, a, 4*4*sizeof(double));
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            if(a[i][j] < num) out[i][j] = num;
}

//FALTA: mvnpdf (Matlab)
/*
    birth_intensity = @(x) (0.1 * mvnpdf(x(1:2), birth_mean1(1:2), covariance_birth(1:2,1:2)) + 0.1 * mvnpdf(x(1:2),birth_mean2(1:2), covariance_birth(1:2,1:2)));//Generate birth weight. This only takes into account the position, not the velocity, as Vo&Ma don't say if they use velocity and I assume that they don't. Taken from page 8 of their paper.
    spawn_intensity = @(x, targetState) 0.05 * mvnpdf(x, targetState, covariance_spawn);//Spawn weight, from page 8 of Vo&Ma.
 */
