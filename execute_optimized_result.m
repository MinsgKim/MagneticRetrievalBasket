%% execute_optimizated_result
clear;clc;
Sim = Call_MagneticBasket_Params;
open("Optimization_250328.slx")

%% Optimization
% Sim.Angle = deg2rad([25 30 45 -45 -35 -25]);
% Sim.Angle = deg2rad([90 90 90 -90 -90 -90]);
% Sim.Angle = deg2rad([84.4715   80.9669   79.3928  -75.4433  -78.0180  -83.8258]);
% Sim.Angle = deg2rad([85 80 80 -75 -80 -85]);
% Sim.Angle = deg2rad(linspace(10, 170, 6));
% Sim.Angle = deg2rad([73.1960   62.4226   57.4241  -56.9934  -62.2620  -73.2091]);
% Sim.Angle = deg2rad([75   60   55  -55  -60  -75]);
Sim.Angle = deg2rad([76.7735   67.8253   62.0523  -61.3219  -67.3133  -76.6115]);

% Sim.Stiffness = 2.0362e-6;
% Sim.Stiffness_upper = 1.8320e-6;

Sim.Plate_MagnetMagnitude = ones(1,6)*120000*Sim.PlateVolume;   % magnetization of each link
Sim.Stiffness = 2.5147e-6;                                      % spring constant
Sim.Stiffness_upper = 2.2831e-6;                                % spring constant for constraints

Sim.BField_Magnitude = 20;      % External B-field with respect to -x axis
Sim.BField_tilt = 0;            % External B-field with respect to +y axis
Sim.Transition_upper = 0.1745;

% Sim.Plate_MagnetMagnitude = [13 13 13 13 13 13] * 1e-4;
% Sim.Stiffness = 6.5131e-6;
% Sim.Stiffness_upper = 9.3752e-5;
% MB.Damping = 1e-7*2; % [Nm/(radian/s)]
% MB.Damping_upper = 2e-7*2; % [Nm/(radian/s)]

[Upper_q, Lower_q] = Sim_Optimization(Sim, "GUI", "on", "DispInfo", "yes");

[Fitness_, v] = optim_bending_and_stability_disp(Sim.Angle, Sim, "GUI", "off", "DispInfo", "no", "InputAngleUnit", "Radian");

disp(v)