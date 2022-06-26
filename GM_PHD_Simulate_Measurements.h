//
// Created by david on 25/05/2022.
//

#ifndef PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H
#define PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H

#define NUM_F 10
#define NUM 30000

void matrizVector(double mat[][4], int n, double vec[], int m, double out[], int &nout);
bool isempty(double v[], int n);
void GM_PHD_Simulate_Measurements(void);


#endif //PROYECTO_GM_PHD_SIMULATE_MEASUREMENTS_H
