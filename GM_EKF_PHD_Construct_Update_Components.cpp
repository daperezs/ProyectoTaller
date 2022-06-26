//
// Created by david
//

#include <iostream>
#include <cstring>
#include "GM_EKF_PHD_Construct_Update_Components.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"

using namespace std;

void cholCUC(double matrix[][4], int n, int m, double out[][4]);
void eyeCUC(double out[][4]);
void matrizxMatrizCUC(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc);

void GM_EKF_PHD_Construct_Update_Components(void)
{
    //Previous define variables
    extern double mk_k_minus_1[NUM_F][NUM], H2[4][4], Pk_k_minus_1[NUM_F][NUM], R2[4][4], x_sensor[3], eta[NUM], S[4][NUM],
            K[4][NUM], P_k_k[4][NUM], m_j[4], eta_j[4], PHt[4][4], S_j[4][4], SChol[4][4], SCholInv[4][4],
            K_j[4][4], W1[4][4], P_j[4][4], sigma_r, sigma_theta, H[2][4], dt, U[4][4];;
    extern int numTargets_Jk_k_minus_1, P_range[4], ceta, cS, cK, cP_k_k;

    //New variables
    extern double xR, yR, hR, R_polar[2][2], r, theta, R_cartesian[2][2], R_velocity[2][2], R[4][4], HPHt[4][4];

    int n, nf, nc;
    double temp[4], out[4], temp1[4][4], temp2[2][2], temp3[4][4];


    //This file creates the components needed for performing a Kalman filter update on the
    //targets using the measurement.
    cout << "Step 3: Constructing update components for all targets, new and existing." << endl;

    //We need to clear the data structures each iteration

    ceta = 0;
    cS = 0;
    cK = 0;
    cP_k_k = 0;

    xR = x_sensor[1];
    yR = x_sensor[2];
    hR = x_sensor[3];
    calculate_R_polar(sigma_r, sigma_theta, R_polar, n);//Sensor covariance, polar


    for(int j = 0; j< numTargets_Jk_k_minus_1; j++)
    {
        for(int i = 0; i<4; i++)
            m_j[i] = mk_k_minus_1[i][j];

        memset(eta_j,0, sizeof(eta_j));//Expected observation
        h(xR,yR,hR,m_j[1], m_j[2], eta_j, n);//Observation model.
        eta_j[2] = m_j[2];
        eta_j[3] = m_j[3]; //Assume constant velocity.
        calculateDataRange4(j, P_range, n); //4x4 array

        calculate_Jacobian_H(xR, yR, m_j[1], m_j[2], H, n);

        for(int i=0; i<4; i++)
            for(int j=0; j<4; j++)
                temp1[i][j] = Pk_k_minus_1[i][j];

        matrizxMatrizCUC(temp1,4,4,H,4,4,PHt, nf,nc); //Taken from Tim Bailey's EKF code. 4x2 array

        r = eta_j[1];
        theta = eta_j[2];
        //We need to extend the sensor covariance R to include the covariance of
        //the speed estimation. See GM_PHD_Initialise_Jacobians for details on
        //this calculation.
        calculate_R_cartesian(R_polar[0][0], r, theta, hR, R_cartesian, n); //sensor covariance, cartesian

        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
                temp2[i][j] = Pk_k_minus_1[i][j];

        //calculate_R_velocity_cartesian(temp2, 2, R_cartesian, 2, dt, R_velocity, n); //We need the landmark position uncertainty, the sensor uncertainty (in cartesian space) and the time delta

        memset(R_polar,0, 4*4*sizeof(double));
        R[0][0] = R_polar[0][0];
        R[0][1] = R_polar[0][1];
        R[1][0] = R_polar[1][0];
        R[1][1] = R_polar[1][1];

        R[2][3] = R_velocity[0][0];
        R[2][3] = R_velocity[0][1];
        R[3][3] = R_velocity[1][0];
        R[3][3] = R_velocity[1][1];

        //Calculate K via Tim Bailey's method.
        //U is set in GM_PHD_Initialise_Jacobians.
        matrizxMatrizCUC(U,4,4,R,4,4,temp1, nf,nc);
        matrizxMatrizCUC(temp1,4,4,U,4,4,temp1, nf,nc);
        matrizxMatrizCUC(H,4,4,PHt,4,4,temp3, nf,nc);
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
                S_j[i][j] = temp1[i][j] + temp3[i][j];

        //At this point, Tim Bailey's code makes S_j symmetric. In this case, it leads to the matrix being non-positive definite a lot of the time and chol() crashes.
        //So we won't do that.
        cholCUC(S_j,4,4, SChol);

        eyeCUC(SChol);
        for(int k=0; k<4; k++)
            SCholInv[k][j] = SChol[k][j] / SChol[k][j]; // triangular matrix, invert via left division

        matrizxMatrizCUC(PHt, 4, 4, SCholInv, 4, 4, W1, nf, nc);

        matrizxMatrizCUC(H,4,4,PHt,4,4,temp1, nf,nc);
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
                HPHt[i][j] = temp1[i][j];

        matrizxMatrizCUC(W1, 4, 4, SCholInv, 4, 4, K_j, nf, nc);

        for(int k=0; k<4; k++)
        {
            temp[k] = Pk_k_minus_1[k][P_range[j]];
        }
        matrizVector(H2, 4, temp, 4, out, n);
        matrizxMatrizCUC(W1, 4, 4, W1, 4, 4, W1, nf, nc);
        for(int k=0; k<4; k++)
        {
            P_j[k][j] = R2[k][j] - W1[k][j];//4x4 array
        }
        //End Tim Bailey's code.

        for(int k = 0; k<4; k++)
        {
            eta[ceta + k] = eta_j[k];

            S[0][cS + k] = S_j[0][k];
            S[1][cS + k] = S_j[1][k];
            S[2][cS + k] = S_j[2][k];
            S[3][cS + k] = S_j[3][k];

            K[0][cK + k] = K_j[0][k];
            K[1][cK + k] = K_j[1][k];
            K[2][cK + k] = K_j[2][k];
            K[3][cK + k] = K_j[3][k];

            P_k_k[0][cP_k_k + k] = P_j[0][k];
            P_k_k[1][cP_k_k + k] = P_j[1][k];
            P_k_k[2][cP_k_k + k] = P_j[2][k];
            P_k_k[3][cP_k_k + k] = P_j[3][k];

            ceta = ceta + k + 1;
            cS = cS + k + 1;
            cK = cK + k + 1;
            cP_k_k = cP_k_k + k + 1;
        }
    }
}



void cholCUC(double matrix[][4], int n, int m, double out[][4])
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (i > j)
            {
                out[i][j] = 0;
            }
            else
                out[i][j] = matrix[i][j];
        }
        cout << endl;
    }
}

void eyeCUC(double out[][4])
{
    for (int c = 0; c < 4; c++) {
        for (int f = 0; f < 4; f++) {
            if (f == c)
                out[f][c] = 1;
            else
                out[f][c] = 0;
        }
    }
}

void matrizxMatrizCUC(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc)
{
    double sum=0;
    int i, j, k;

    if(ac != bf)
    {
        std::cout << "Different dimensions\n";
        exit(0);
    }
    nf = af;
    nc = bc;

    for(i=0;i<af;i++)
    {
        for(j=0;j<bf;j++)
        {
            for(k=0;k<bc;k++)
            {
                sum+=m1[i][k]*m2[k][j];
            }
            out[i][j]=sum;
            sum=0;
        }
    }
}