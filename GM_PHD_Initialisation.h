//
// Created by david on 04/05/2022.
//

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
