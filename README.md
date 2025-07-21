# Satellite Attitude Control

This repository contains a complete MATLAB/Simulink project focused on the modeling, simulation, and control of the attitude dynamics of a 3U CubeSat. The work is structured into three progressive phases aimed at developing, tuning, and analyzing a full Attitude Determination and Control System (ADCS), incorporating both ideal and non-ideal effects.

## Objectives

1. **Attitude Dynamics and Kinematics Simulation**  
   - Development of a Simulink model implementing the translational and rotational dynamics of a CubeSat.
   - Input: sinusoidal torque profiles.  
   - Output: satellite attitude expressed in both quaternion representation and rotation matrix.

2. **Controller Design and Integration**  
   - Design of a feedback controller using the *MATLAB Control System Designer* app.
   - Controller parameters are dimensioned based on the satellite's physical properties and implemented via a separate MATLAB script.
   - Integration of the controller with the Simulink model from Phase 1.

3. **Realistic CubeSat Attitude Control Simulation**  
   - Implementation of a full 3-axis stabilized, nadir-pointing CubeSat model including:
     - External disturbances: gravity gradient, Earth’s magnetic field, etc.
     - Actuation systems: reaction wheels and magnetorquers with non-ideal behavior.
     - Sensor models: noisy measurements and limited accuracy.
     - Solar panels and structural effects.
   - Disturbance rejection techniques for low-frequency perturbations.
   - Frequency-domain analysis using Fourier Transform of the attitude response.

