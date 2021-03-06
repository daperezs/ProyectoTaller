//
// Created by david
//

/*
%GM_EKF_PHD_Create_Birth
%Matlab code by Bryan Clarke b.clarke@acfr.usyd.edu.au

%This file performs the processing on measurements to extract new targets
%and populates the birth lists (mean, weight, covariance) needed to instantiate them next iteration.
*/

#include <iostream>
#include <cstring>
#include "GM_EKF_PHD_Create_Birth.h"

using namespace std;

void  GM_EKF_PHD_Create_Birth(void)
{
    //Previos define variables
    extern double w_birth[NUM],m_birth[NUM_F][NUM], P_birth[NUM_F][NUM], w_spawn[NUM], m_spawn[NUM_F][NUM],
            P_spawn[NUM_F][NUM], k, thisMeasRowRange, prevMeasRowRange;;
    extern int numBirthedTargets, numSpawnedTargets, addVelocityForNewTargets;


    //This is not formally spelt out in Vo&Ma so I have used my best judgement to
    //come up with a method to initialise targets. Targets are initialised as
    //static (using the position of current measurements) and/or dynamic (using
    //two generations of measurements and calculating the movement between).
    //If only static targets are added, the number of new targets = the number
    //of measurements.
    //If only dynamic targets are added, the number of new targets = the number of
    //measurements this teration * the number of measurements last iteration
    //If both static and dynamic targets are added, the number of new targets is
    //equal to the sum of these.

    cout << "Step 7: Creating new targets from measurements, for birthing next iteration" << endl;
    memset(w_birth,0, sizeof(w_birth));
    memset(m_birth,0, NUM_F*NUM*sizeof(double));
    memset(P_birth,0, NUM_F*NUM*sizeof(double));
    memset(w_spawn,0, sizeof(w_spawn));
    memset(m_spawn,0, NUM_F*NUM*sizeof(double));
    memset(P_spawn,0, NUM_F*NUM*sizeof(double));
    numBirthedTargets = 0;
    numSpawnedTargets = 0;
}

/*



%We create targets using two generations of measurements.
%The first generation is used to calculate the velocity (dx/dt) to the
%second generation.
%The second generation gives the position.
%We also add a set of targets from the second generation but with velocity
%zero, since some may be unlinked to the first generation, but we have no
%velocity information for them.
if((addVelocityForNewTargets == true) && (k >= 2))%If we want to add targets with initial velocities.If only one iteration complete, cannot calculate velocity
    %Each measurement consists of 2 rows
    thisMeasRowRange = k;
    prevMeasRowRange = k-1;

    thisMeas = simMeasurementHistory{thisMeasRowRange};
    prevMeas = simMeasurementHistory{prevMeasRowRange};

    for j_this = 1:size(thisMeas,2)
        for j_prev = 1:1:size(prevMeas,2)%Match every pair from previous to current
            z_this = thisMeas(:,j_this);
            z_prev = prevMeas(:, j_prev);
            %Calculate and add the velocity.
            m_i = inv_h(z_this(1), z_this(2), x_sensor(1), x_sensor(2), x_sensor(3));
            thisV = (m_i(1:2) - inv_h(z_prev(1), z_prev(2), x_sensor(1), x_sensor(2), x_sensor(3))) / dt;
            if(abs(thisV(1)) > MAX_V) || (abs(thisV(2)) > MAX_V)
                continue;%To reduce the number of targets added, we filter out the targets with unacceptable velocities.
            end
            %Augment with velocity
            m_i = [m_i; thisV];

            %Decide if the target is birthed (from birth position)
            %or spawned (from an existing target)
            %Initialise the weight to birth
            birthWeight = birth_intensity(m_i);
            %Targets can also spawn from existing targets. We will
            %take whichever is a higher weight - birthing or
            %spawning
            nTargets = size(X_k, 2);
            maxSpawnWeight = -1;
            for targetI = 1:nTargets
                thisWeight = spawn_intensity(m_i, X_k(:,targetI)) * X_k_w(targetI);%Spawn weight is a function of proximity to the existing target, and the weight of the existing target.
                if(thisWeight > maxSpawnWeight)
                    maxSpawnWeight = thisWeight;
                end
            end
            %Check if birthing had higher weight.
            if(birthWeight > maxSpawnWeight)
                %Birth the target
                w_i = birthWeight;
                %Initialise the covariance
                P_i = covariance_birth;
                w_birth = [w_birth, w_i];
                m_birth = [m_birth m_i];
                P_birth = [P_birth, P_i];
                numBirthedTargets = numBirthedTargets + 1;
            else
                %Spawn the target
                w_i = maxSpawnWeight;
                %Initialise the covariance
                P_i = covariance_spawn;
                w_spawn = [w_spawn, w_i];
                m_spawn = [m_spawn, m_i];
                P_spawn = [P_spawn, P_i];
                numSpawnedTargets = numSpawnedTargets + 1;
            end
        end
    end
end
%If we want to add targets, treating them as if they are
%static.
if (addStaticNewTargets == true)
     thisMeasRowRange = k;
    thisMeas = simMeasurementHistory{thisMeasRowRange};
    for j_this = 1:size(thisMeas,2)    %Each measurement consists of 2 rows

        %Add a static target
        z_this = thisMeas(:,j_this);
        m_i = inv_h(z_this(1), z_this(2), x_sensor(1), x_sensor(2), x_sensor(3));

        m_i(3:4) = [0; 0];

        %Decide if the target is birthed (from birth position)
        %or spawned (from an existing target)
        %Initialise the weight to birth
        birthWeight = birth_intensity(m_i);
        %Targets can also spawn from existing targets. We will
        %take whichever is a higher weight - birthing or
        %spawning
        nTargets = size(X_k, 2);
        maxSpawnWeight = -1;
        for targetI = 1:nTargets
            thisWeight = spawn_intensity(m_i, X_k(:,targetI)) * X_k_w(targetI);%Spawn weight is a function of proximity to the existing target, and the weight of the existing target.
            if(thisWeight > maxSpawnWeight)
                maxSpawnWeight = thisWeight;
            end
        end
        %Check if birthing had higher weight.
        if(birthWeight > maxSpawnWeight)
            %Birth the target
            w_i = birthWeight;
            %Initialise the covariance
            P_i = covariance_birth;
            w_birth = [w_birth, w_i];
            m_birth = [m_birth m_i];
            P_birth = [P_birth, P_i];
            numBirthedTargets = numBirthedTargets + 1;
        else
            %Spawn the target
            w_i = maxSpawnWeight;
            %Initialise the covariance
            P_i = covariance_spawn;
            w_spawn = [w_spawn, w_i];
            m_spawn = [m_spawn m_i];
            P_spawn= [P_spawn, P_i];
            numSpawnedTargets = numSpawnedTargets + 1;
        end
    end
end


if VERBOSE == 1
    for j = 1:numBirthedTargets
        thisM = m_birth(:,j);
        s = sprintf('Target to birth %d: %3.4f %3.4f %3.4f %3.4f Weight %3.9f', j, thisM(1), thisM(2), thisM(3), thisM(4), w_birth(j));
        disp(s);
    end
    for j = 1:numSpawnedTargets
        thisM = m_spawn(:,j);
        s = sprintf('Target to spawn %d: %3.4f %3.4f %3.4f %3.4f Weight %3.9f', j, thisM(1), thisM(2), thisM(3), thisM(4), w_spawn(j));
        disp(s);
    end
end

 */