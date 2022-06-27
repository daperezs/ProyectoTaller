//
// Created by david
//

/**
    @file GM_EKF_PHD_Simulate_Measurements.h
    @title GM_EKF_PHD_Simulate_Measurements
    @brief This file generates simulated measurement data for the linear KF simulation
    described in Vo&Ma. There is a nonlinear simulation described in the paper
    but I have not gotten around to implementing it.
    The sensor measurements are range-bearing measurements, similar to those
    from an FMCW radar. The expected error magnitude is independent of
    measurement (i.e. range error is not proportional to range, bearing error
    is not proportional to bearing).
 */

#ifndef PROYECTO_GM_EKF_PHD_SIMULATE_MEASUREMENTS_H
#define PROYECTO_GM_EKF_PHD_SIMULATE_MEASUREMENTS_H

#define NUM_F 10
#define NUM 30000

void GM_EKF_PHD_Simulate_Measurements(void);


#endif //PROYECTO_GM_EKF_PHD_SIMULATE_MEASUREMENTS_H
