//
// Created by david
//

/**
    @file GM_EKF_PHD_Update.h
    @title GM_EKF_PHD_Update
    @brief This file performs a PHD filter update on the targets. This is basically a brute-force Kalman update of every target with every measurement and creating a new target from the update results.
*/

#ifndef PROYECTO_GM_EKF_PHD_UPDATE_H
#define PROYECTO_GM_EKF_PHD_UPDATE_H

#define NUM_F 10
#define NUM 30000

void GM_EKF_PHD_Update(void);


#endif //PROYECTO_GM_EKF_PHD_UPDATE_H
