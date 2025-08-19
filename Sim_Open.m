%%
clear; clc;
Params = Sim_Param_Setup;

ParamsBus = Simulink.Bus.createObject(Params);

open('tests.slx')

% tests.slx: one-arm simulation
% Copy_of_tests.slx: mult-arm simulation (progressing)