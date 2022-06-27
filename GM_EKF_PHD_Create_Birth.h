//
// Created by david
//

/**
    @file GM_EKF_PHD_Create_Birth.h
    @title GM_EKF_PHD_Create_Birth
    @brief This is not formally spelt out in Vo&Ma so I have used my best judgement to
    come up with a method to initialise targets. Targets are initialised as
    static (using the position of current measurements) and/or dynamic (using
    two generations of measurements and calculating the movement between).
 */

#ifndef PROYECTO_GM_EKF_PHD_CREATE_BIRTH_H
#define PROYECTO_GM_EKF_PHD_CREATE_BIRTH_H

#define NUM_F 10
#define NUM 30000

void  GM_EKF_PHD_Create_Birth(void);


#endif //PROYECTO_GM_EKF_PHD_CREATE_BIRTH_H
