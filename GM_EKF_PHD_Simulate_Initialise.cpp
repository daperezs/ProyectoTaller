//
// Created by david
//

//GM_EKF_PHD_Simulate_Initialise
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

//We implement the simulation described in example 1 of Vo&Ma 2006, the one
//designed for the linear KF PHD filter. Although a second, different
//example is listed for the EKF version, I haven't gotten around to
//implementing it yet.

//Most of this is unchanged from the example used in the linear KF PHD
//filter. The big changes for the EKF are in
//GM_EKF_PHD_Initialise_Jacobians.

#include <cstring>
#include "GM_EKF_PHD_Simulate_Initialise.h"

void GM_EKF_PHD_Simulate_Initialise(void){

    // Previous define
    extern int nClutter, endTime, c1History, c2History, c3History;
    extern double noiseScaler, simTarget1Start[4], simTarget2Start[4], simTarget3Start[4],simTarget1End[2],
            simTarget2End[2], simTarget3End[2], simTarget1Vel[2], simTarget2Vel[2], simTarget3Vel[2],
            simTarget1History[NUM_F][NUM], simTarget2History[NUM_F][NUM], simTarget3History[NUM_F][NUM],
            simTarget1State[4], simTarget2State[4], simTarget3State[4], simTarget3SpawnTime, birth_mean1[4],
            birth_mean2[4];

    // New global variables
    extern double x_sensor[3];

    // Control parameters
    memset(x_sensor,0, sizeof(x_sensor)); //The position of the simulated range-bearing sensor. Heading 0 is along the positive X-axis. Counter clockwise is positive, clockwise is negative, bearings are bounded -pi < theta <= pi
    noiseScaler = 1.0; //Adjust the strength of the noise on the measurements by adjusting this
    nClutter = 50; //Assume constant 50 clutter measurements. Since clutter is Poisson distrbuted it might be more accurate to use nClutter = poissrnd(50) if you have the required Matlab toolbox. Constant 50 clutter works well enough for simulation purposes.

    // Target parameters
    endTime = 100; //Duration of main loop.
    //simTarget1Start = birth_mean1;
    memcpy(simTarget1Start, birth_mean1, 4*sizeof(double));
    //simTarget2Start = birth_mean2;
    memcpy(simTarget2Start, birth_mean2, 4*sizeof(double));

    simTarget1End[0] = 500;
    simTarget1End[1]= -900;
    simTarget2End[0] = 900;
    simTarget2End[1]= -500;
    simTarget3End[0] = -200;
    simTarget3End[1]= -750;

    simTarget1Vel[0] = (simTarget1End[0] - simTarget1Start[0]) / endTime;
    simTarget1Vel[1] = (simTarget1End[1] - simTarget1Start[1]) / endTime;
    simTarget2Vel[0] = (simTarget2End[0] - simTarget2Start[0]) / endTime;
    simTarget2Vel[1] = (simTarget2End[1] - simTarget2Start[1]) / endTime;
    simTarget3Vel[0] = -20;
    simTarget3Vel[1]= -4;

    simTarget1Start[2] = simTarget1Vel[0];
    simTarget1Start[3] = simTarget1Vel[1];
    simTarget2Start[2] = simTarget2Vel[0];
    simTarget2Start[3] = simTarget2Vel[1];
    simTarget3Start[0] = simTarget3Start[1] = 0;
    simTarget3Start[2] = simTarget3Vel[0];
    simTarget3Start[3] = simTarget3Vel[1];

    c1History = 0;
    c2History = 0;
    c3History = 0;
    memset(simTarget1History,0, NUM_F*NUM*sizeof(double));
    memset(simTarget2History,0, NUM_F*NUM*sizeof(double ));
    for(int i=0; i<4; i++)
    {
        simTarget1History[i][0] = simTarget1Start[i];
        simTarget2History[i][0] = simTarget2Start[i];
    }
    c1History++;
    c2History++;
    memset(simTarget3History,0, 4*NUM*sizeof(double));

    /* FALTA DEFINIR
    simTrueMeasurementHistory = [];
    simNoisyMeasurementHistory = {}; //We use a cell array so that we can have rows of varying length.
    simVehicleStateHistory = {}; //We store a history of vehicle positions at the same time as measurements are made.
     */

    //simTarget2Start = simTarget1Start;
    memcpy(simTarget1State, simTarget1Start, 4*sizeof(double));
    //simTarget2State = simTarget2Start;
    memcpy(simTarget2State, simTarget2Start, 4*sizeof(double));
    //simTarget3State = [];
    memset(simTarget3State,0, sizeof(simTarget3State));

    simTarget3SpawnTime = 66; //Target 3 is spawned from target 1 at t = 66s
}
/* FALTA

%Set up for plot
%Measurements and targets plot
figure(1);
clf;
hold on;
axis([-1000 1000 -1000 1000]);
xlim([-1000 1000]);
ylim([-1000 1000]);

%X and Y measurements plot
figure(2);
clf;
subplot(2,1,1);
hold on;
axis([0, 100, -1000, 1000]);
xlabel('Simulation step');
ylabel('X position of measurement (m)');
title('Measurement X coordinates');
subplot(2,1,2);
hold on;
axis([0, 100, -1000, 1000]);
xlabel('Simulation step');
ylabel('Y position of measurement (m)');
title('Measurement Y coordinates');

%Performance metric plot
if(CALCULATE_OSPA_METRIC == 1)
    figure(3);
    clf;
    xlabel('Simulation step');
    ylabel('OSPA error metric (higher is worse)');
end
 */

