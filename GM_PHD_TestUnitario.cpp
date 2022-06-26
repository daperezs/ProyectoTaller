//
// Created by david on 05/05/2022.
//

#include <iostream>
#include <stdio.h>
#include <cmath>
#include "globales.h"
#include "GM_funcionAnonima.h"
#include "unifpdf_2d.h"
#include "GM_EKF_PHD_Initialise_Jacobians.h"
#include "GM_PHD_Simulate_Measurements.h"

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int tests_run = 0;

int test();
int all_tests();
bool comparaV(int a[], int na, int b[], int nb);
bool comparaV(double a[], int na, double b[], int nb);
bool comparaM(double a[][4], int na, double b[][4], int nb);
bool comparaM(double a[][2], int na, double b[][2], int nb);
int calculateDataRange2_01();
int calculateDataRange4_01();
int unifpdf_2d_01();
int clutter_intensity_01();
int hypot_01();
int h_01();
int inv_h_01();
int calculate_R_polar_01();
int Jacobian_of_velocity_f_01();
int calculate_Jacobian_of_inv_h_01();
int traspuesta_01();
int matrizMatriz_01();
int calculate_R_cartesian_01();
int calculate_R_velocity_cartesian_01();
int calculate_Jacobian_F_01();
int calculate_Jacobian_G_01();
int matrizVector_01();
int isemptyEKF_01();
int chol_01();

int testUnitario()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    cout << "Tests run: " << tests_run << endl;

    return result != 0;
}

int all_tests()
{
    _verify(calculateDataRange2_01);
    _verify(calculateDataRange4_01);
    _verify(unifpdf_2d_01);
    _verify(clutter_intensity_01);
    _verify(hypot_01);
    _verify(h_01);
    _verify(inv_h_01);
    _verify(calculate_R_polar_01);
    _verify(Jacobian_of_velocity_f_01);
    _verify(calculate_Jacobian_of_inv_h_01);
    _verify(traspuesta_01);
    _verify(matrizMatriz_01);
    _verify(calculate_R_cartesian_01);
    _verify(calculate_R_velocity_cartesian_01);
    _verify(calculate_Jacobian_F_01);
    _verify(calculate_Jacobian_G_01);
    _verify(matrizVector_01);
    _verify(isemptyEKF_01);
    _verify(chol_01);

    return 0;
}

bool comparaV(int a[], int na, int b[], int nb)
{
    if(na != nb) return false;
    for(int i=0; i<na; i++)
        if(a[i] != b[i]) return false;

    return true;
}

bool comparaV(double a[], int na, double b[], int nb)
{
    if(na != nb) return false;
    for(int i=0; i<na; i++)
        if(a[i] != b[i]) return false;

    return true;
}

bool comparaM(double a[][4], int na, double b[][4], int nb)
{
    int k=0;

    if(na != nb) return false;

    for(int i=0; i<na; i++)
        for(int j=0; j<nb; j++)
            if(a[i][j] != b[i][j]) return false;

    return true;
}

int calculateDataRange2_01()
{
    int out[2], n, test[2]={111,112};

    calculateDataRange2(56, out, n);

    cout << "Test calculateDataRange2" << endl;
    cout << out[0] << endl;
    cout << out[1] << endl;

    _assert(comparaV(out, n, test, n));

    return 0;
}

int calculateDataRange4_01()
{
    int out[4], n, test[4]={309,310,311,312};

    calculateDataRange4(78, out, n);

    cout << "Test calculateDataRange2" << endl;
    cout << out[0] << endl;
    cout << out[1] << endl;
    cout << out[2] << endl;
    cout << out[3] << endl;

    _assert(comparaV(out, n, test, n));

    return 0;
}

int unifpdf_2d_01()
{
    int xrange[] = {-1000,1000}, yrange[] = {-1000,1000};
    double z[]= {100,100}, val = 2.5000e-07;

    cout << "Test unifpdf_2d" << endl;

    _assert(unifpdf_2d(xrange, yrange, z) - val < pow(10,-12));

    return 0;
}

int clutter_intensity_01()
{
    int xrange[] = {-1000,1000}, yrange[] = {-1000,1000};
    double z[]= {100,100}, val = 2.5000e-05;

    cout << "Test clutter_intensity" << endl;

    _assert(clutter_intensity(z) - val < pow(10,-12));

    return 0;
}

int hypot_01()
{
    double val = 1.4142;

    cout << "Test hypot" << endl;

    _assert(hypot(1,12) - val < pow(10,-12));

    return 0;
}

int h_01()
{
    double out[2], test[] = {0,-1};
    int n;

    cout << "Test hypot" << endl;

    h(1.0,1.0,1.0,1.0,1.0, out, n);

    _assert(comparaV(out, n, test, n));

    return 0;
}

int calculate_Jacobian_H()
{
    double out[4][4], test[4][4]= {{1.0,0,0,0},
                                            {0,1.0,0,0},
                                            {0,0,1.0,0},
                                            {0,0,0,1.0}};
    int n;

    cout << "Test calculate_Jacobian" << endl;

    calculate_Jacobian_H(1.0,1.0,1.0,1.0, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int inv_h_01()
{

    double out[2], test[]= {0.5839,1.9093};
    int n;

    cout << "Test inv_h" << endl;

    inv_h(1.0,1.0,1.0,1.0,1.0, out, n);

    _assert(comparaV(out, n, test, n));

    return 0;
}

int calculate_R_polar_01()
{

    double out[2][2], test[][2]= {{1.0,0},
                                   {0,1.0}};
    int n;

    cout << "Test calculate_R_polar" << endl;

    calculate_R_polar(1.0,1.0, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int Jacobian_of_velocity_f_01()
{
    double out[4][4], test[4][4]= {{-1.0,1.0,0,0},
                                   {0,0,-1.0,1.0}};
    int n, m;

    cout << "Test Jacobian_of_velocity_f" << endl;

    Jacobian_of_velocity_f(1.0, out, n, m);

    _assert(comparaM(out, n, test, n));

    return 0;

}

int calculate_Jacobian_of_inv_h_01()
{
    double out[2][2], test[][2]= {{-0.4161,-0.9093},
                                  {0.9093,-0.4161}};
    int n;

    cout << "Test calculate_Jacobian_of_inv_h" << endl;

    calculate_Jacobian_of_inv_h(1.0,1.0,1.0, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;

}

int traspuesta_01()
{
    double out[2][2], test[][2]= {{1,1},
                                  {2,2}};
    double mat[][2] = {{1,2},
                       {1,2}};
    int n;

    cout << "Test traspuesta" << endl;

    traspuesta(mat, 2, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;

}

int matrizMatriz_01()
{
    double out[2][2], test[][2]= {{3,6},
                                  {3,6}};
    double m1[][2] = {{1,2},
                       {1,2}};
    double m2[][2] = {{1,2},
                      {1,2}};
    int n , m;

    cout << "Test matrizMatriz" << endl;

    matrizMatriz(m1, 2,2, m2, 2,2, out, n, m);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int calculate_R_cartesian_01()
{
    double out[2][2], test[][2]= {{1,0},
                                  {0,1}};
    int n , m;

    cout << "Test calculate_R_cartesian" << endl;

    calculate_R_cartesian(1,1,1,1, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int calculate_R_velocity_cartesian_01()
{
    void calculate_R_velocity_cartesian(double P_landmark[][2], int n, double R_c[][2], int m, double dt, double outt[][2], int no);
    double out[2][2], test[][2]= {{1,0},
                                  {0,1}};
    int n;

    double m1[][2] = {{1,2},
                      {1,2}};
    double m2[][2] = {{1,2},
                      {1,2}};

    cout << "Test calculate_R_velocity_cartesian" << endl;

    calculate_R_velocity_cartesian(m1,1,m2,1,1, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int calculate_Jacobian_F_01()
{
    double out[4][4], test[4][4]= {{1.0,0,1.0,0},
                                   {0,1.0,0,1.0},
                                   {0,0,1.0,0},
                                   {0,0,0,1.0}};
    int n;

    cout << "Test calculate_Jacobian_F" << endl;

    calculate_Jacobian_F(1.0, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int calculate_Jacobian_G_01()
{
    double out[4][4], test[4][4]= {{1.0,0,1.0,0},
                                   {0,1.0,0,1.0},
                                   {0,0,1.0,0},
                                   {0,0,0,1.0}};
    int n;

    cout << "Test calculate_Jacobian_G" << endl;

    calculate_Jacobian_G(1.0, out, n);

    _assert(comparaM(out, n, test, n));

    return 0;
}

int matrizVector_01()
{

    double out[4], test[4]= {2,2,1,1};
    int n;

    double m[4][4]= {{1.0,0,1.0,0},
                     {0,1.0,0,1.0},
                     {0,0,1.0,0},
                     {0,0,0,1.0}};
    double v[4]={1,1,1,1};

    cout << "Test matrizVector" << endl;

    matrizVector(m, 4, v, 4,out, n);

    _assert(comparaV(out, n, test, n));

    return 0;

}

int isemptyEKF_01()
{
    bool isemptyEKF(double v[], int n);

    double out[4];
    int n;


    cout << "Test matrizVector" << endl;

    bool aux = isempty(out, n);

    _assert(aux = true);

    return 0;
}

int chol_01()
{
    void chol(double matrix[][4], int n, int m, double out[][4]);

    double out[4][4], test[4][4]= {{1.0,0,1.0,0},
                                   {0,1.0,0,1.0},
                                   {0,0,1.0,0},
                                   {0,0,0,1.0}};
    int n;
    double m[4][4]= {{1.0,0,1.0,0},
                     {0,1.0,0,1.0},
                     {0,0,1.0,0},
                     {0,0,0,1.0}};

    cout << "Test calculate_Jacobian_G" << endl;

    chol(m,1,1, out);

    _assert(comparaM(out, n, test, n));

    return 0;
}
