%% execute_calibrated_result
clear;clc;
Sim = Call_MagneticBasket_Params;
open("Calibration_250328.slx")

%% Calibration
% Calibration은 어떻게 실험할건지 정해서 그에 맞게 모델을 변경해야함
% Sim.Plate_MagnetMagnitude = [25004 25001 26963 25000 26339
% 28469]*Sim.PlateVolume; % 1st simulation results Kt_low = 2e-6, Kt_up =
% 3e-6

% Sim.Plate_MagnetMagnitude = [31105 30247 39991 30001 31387 33941]*Sim.PlateVolume;
% Sim.Plate_MagnetMagnitude = [25031 30227 39991 30075 31391 33911]*Sim.PlateVolume;
% Sim.Stiffness = 2.0362e-6;
% Sim.Stiffness_upper = 1.8320e-6;

Sim.Plate_MagnetMagnitude = [30973 30962 39510 25856 25211 25256]*Sim.PlateVolume;
Sim.Stiffness = 2.5147e-6;
Sim.Stiffness_upper = 2.2831e-6;

% Sim.Plate_MagnetMagnitude = [34610 33890 39270 28910 39680 31830]*Sim.PlateVolume;
% Sim.Stiffness = 3.9973e-6;
% Sim.Stiffness_upper = 5.9993e-6;

Sim.Source_Angle = -90;
Sim.Source_Distance_x = 50;
Sim.Source_Distance_z = 15;

Sim.Transition_upper = 0.1745; % 10 degree
% Sim.Transition_upper = 0.6; % ?? degree

% Sim.Plate_MagnetMagnitude = [7 7 7 7 13 13] * 1e-4;
% Sim.Stiffness = 6.5131e-6;
% Sim.Stiffness_upper = 9.3752e-5;
% Sim.Damping = 1e-7*2; % [Nm/(radian/s)]
% Sim.Damping_upper = 2e-7*2; % [Nm/(radian/s)]

Sim.Cal_EndTime = 0.5;
Sim.Angle = deg2rad([25 30 45 -45 -35 -25]);
q = Sim_Calibration(Sim, "GUI", "on", "DispInfo", "yes");