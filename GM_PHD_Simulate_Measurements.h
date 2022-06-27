//
// Created by david
//

/**
    @file GM_PHD_Simulate_Measurements.h
    @title GM_PHD_Simulate_Measurements
    @brief This file generates simulated measurement data for  the simulation
*/

#ifndef PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H
#define PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H

#define NUM_F 10
#define NUM 30000

void matrizVector(double mat[][4], int n, double vec[], int m, double out[], int &nout);
bool isempty(double v[], int n);
void GM_PHD_Simulate_Measurements(void);


#endif //PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H
