%---------- Previous Code ----------%
%% Run_SimOpt_Sejun ver1
clear; clc; close;
rmpath magnetic_robot_sim_opt_MK_v1

addpath('MatlabFunctions')
addpath('magnetic_robot_sim_opt_SP_v1')
Simscape_Preprocessing;

%% Test_Run_SimOpt_Sejun ver1
Fitness = Simscape_Fitness(Input.x_init_test2, SimParams, Fixed, "DispInfo", "yes", "InputAngleUnit", "Radian", "GUI", "on");



%% Run_SimOpt_Minseong ver1
clear; clc; close;
rmpath magnetic_robot_sim_opt_SP_v1

addpath('MatlabFunctions')
addpath('magnetic_robot_sim_opt_MK_v1')
MK_Simscape_Preprocessing;

%% Test_Run_Simopt_Minseong ver1
Fitness = MK_Simscape_Fitness(Input.x_init_1, SimParams, Fixed, "DispInfo", "yes", "InputAngleUnit", "Radian", "GUI", "on");