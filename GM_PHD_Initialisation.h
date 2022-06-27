//
// Created by david
//

/**
    @file GM_PHD_Initialisation.h
    @title GM_PHD_Initialisation
    @brief This file initialises most of the variables that we use for the filter.
*/

#ifndef PROYECTO_GM_PHD_INITIALISATION_H
#define PROYECTO_GM_PHD_INITIALISATION_H

#define NUM_F 10
#define NUM 30000

void GM_PHD_Initialisation(void);
void calculateDataRange2(int j, int out[], int &n);
void calculateDataRange4(int j, int out[], int &n);
double clutter_intensity(double z[]);
void max(double a[][4], int n, double num, double out[][4], int &m);


#endif //PROYECTO_GM_PHD_INITIALISATION_H
