//
// Created by david on 19/06/2022.
//

//GM_PHD_Construct_Update_Components
//Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

#include <iostream>
#include "GM_PHD_Construct_Update_Components.h"
#include "GM_PHD_Simulate_Measurements.h"
#include "GM_PHD_Initialisation.h"

using namespace std;

void chol(double matrix[][4], int n, int m, double out[][4]);
void eye(double out[][4]);
void matrizxMatriz(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc);

void GM_PHD_Construct_Update_Components(void)
{

    //Previous define variables
    extern double mk_k_minus_1[NUM_F][NUM], H2[4][4], Pk_k_minus_1[NUM_F][NUM], R2[4][4];
    extern int numTargets_Jk_k_minus_1, P_range[4];

    //New variables
    extern double eta[NUM], S[4][NUM], K[4][NUM], P_k_k[4][NUM], m_j[4], eta_j[4], PHt[4][4], S_j[4][4], SChol[4][4], SCholInv[4][4],
            K_j[4][4], W1[4][4], P_j[4][4];
    extern int ceta, cS, cK, cP_k_k;

    int n, nf, nc;
    double temp[4], out[4];


    //This file creates the components needed for performing a Kalman filter update on the
    //targets using the measurement.
    cout << "Step 3: Constructing update components for all targets, new and existing." << endl;

    //We need to clear the data structures each iteration

    ceta = 0;
    cS = 0;
    cK = 0;
    cP_k_k = 0;

    for(int j = 0; j< numTargets_Jk_k_minus_1; j++)
    {
        for(int i = 0; i<4; i++)
            m_j[i] = mk_k_minus_1[i][j];

        matrizVector(H2, 4, m_j, 4, eta_j, n); //Observation model. Assume we see position AND velocity of the target.

        calculateDataRange4(j, P_range, n); //4x4 array

        //Taken from Tim Bailey's EKF code. 4x4 array
        for(int k=0; k<4; k++)
        {
            temp[k] = Pk_k_minus_1[k][P_range[k]];
        }
        matrizVector(H2, 4, temp, 4, out, n);
        for(int k=0; k<4; k++)
        {
            PHt[k][j] = out[k];
        }

        //Calculate K via Tim Bailey's method.
        for(int k=0; k<4; k++)
        {
            temp[k] = PHt[k][j];
        }
        matrizVector(H2, 4, temp, 4, out, n);
        for(int k=0; k<4; k++)
        {
            S_j[k][j] = R2[k][j] + out[k];
        }
        //At this point, Tim Bailey's code makes S_j symmetric. In this case, it leads to the matrix being non-positive definite a lot of the time and chol() crashes.
        //So we won't do that.
        chol(S_j,4,4, SChol);

        eye(SChol);
        for(int k=0; k<4; k++)
            SCholInv[k][j] = SChol[k][j] / SChol[k][j]; // triangular matrix, invert via left division

        matrizxMatriz(PHt, 4, 4, SCholInv, 4, 4, W1, nf, nc);

        matrizxMatriz(W1, 4, 4, SCholInv, 4, 4, K_j, nf, nc);

        for(int k=0; k<4; k++)
        {
            temp[k] = Pk_k_minus_1[k][P_range[j]];
        }
        matrizVector(H2, 4, temp, 4, out, n);
        matrizxMatriz(W1, 4, 4, W1, 4, 4, W1, nf, nc);
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

void chol(double matrix[][4], int n, int m, double out[][4])
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

void eye(double out[][4])
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

void matrizxMatriz(double m1[][4], int af, int ac, double m2[][4], int bf, int bc, double out[][4], int &nf, int &nc)
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