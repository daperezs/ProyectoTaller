//
// Created by david
//

/**
    @file GM_EKF_PHD_Initialise_Jacobians.h
    @title GM_EKF_PHD_Initialise_Jacobians
    @brief In this file we calculate the Jacobians for the EKF in the PHD filter.
    These are for the motion and observation models, and the noise on each of
    these.
 */

#ifndef PROYECTO_GM_EKF_PHD_INITIALISE_JACOBIANS_H
#define PROYECTO_GM_EKF_PHD_INITIALISE_JACOBIANS_H


void GM_EKF_PHD_Initialise_Jacobians(void);
double hypot(double x, double y);
void h(double xS, double yS, double hS, double xL, double yL, double out[], int &n);
void calculate_Jacobian_H(double xR, double yR, double xL, double yL, double out[][4], int &n);
void inv_h(double r, double theta, double xS, double yS, double hS, double out[], int &n);
void calculate_R_polar(double sigma_r, double sigma_theta, double out[][2], int &n);
void Jacobian_of_velocity_f(double dt, double out[][4], int &n, int &m);
void calculate_Jacobian_of_inv_h(double r, double theta, double rH, double out[][2], int &n);
void traspuesta(double mat[][2], int n, double tras[][2], int &m);
void matrizMatriz(double m1[][2], int af, int ac, double m2[][2], int bf, int bc, double out[][2], int &nf, int &nc);
void calculate_R_cartesian(double R, double r, double theta, double rH, double out[][2], int &n);
void calculate_Jacobian_F(double dt, double out[][4], int &n);
void calculate_Jacobian_G(double dt, double out[][4], int &n);



#endif //PROYECTO_GM_EKF_PHD_INITIALISE_JACOBIANS_H
