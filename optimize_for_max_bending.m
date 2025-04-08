%% Optimization for Max. Bending
clear;clc;
Sim = Call_MagneticBasket_Params;
open("Optimization_250328.slx")

select_stiff_ver = input('remanence or moment?: ', 's');

if strcmpi(select_stiff_ver, 'remanence')

    Sim.Plate_MagnetMagnitude = ones(1,6) * 120000 * Sim.PlateVolume;
    Sim.Stiffness = 2.5147e-6;
    Sim.Stiffness_upper = 2.2831e-6;
    Sim.Damping = 10e-9*2; % [Nm/(radian/s)]
    Sim.Damping_upper = 10e-9*2;
    Sim.Transition_upper = 0.1745;

else

    Sim.Plate_MagnetMagnitude = ones(1,6) * 15 *1e-4;
    Sim.Stiffness = 6.5131e-6;
    Sim.Stiffness_upper = 9.3752e-5;
    Sim.Damping = 1e-7*2; % [Nm/(radian/s)]
    Sim.Damping_upper = 2e-7*2; % [Nm/(radian/s)]

end

x_init = deg2rad([82.1934   75.7876   74.6678  -74.0767  -75.3764  -82.1305]);

Sim.BField_Magnitude = 20;
Sim.BField_tilt = 0;

lb = repmat(-pi,1,6);
ub = repmat(pi,1,6);

options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-8, 'FunctionTolerance', 1e-8, ...
    'MaxFunctionEvaluations', 1e5, 'OptimalityTolerance', 1e-8, ...
    'Algorithm', 'interior-point', "EnableFeasibilityMode", true, ...
    "SubproblemAlgorithm","cg");

tic;
[x_opt, f_val] = fmincon(@(x) optim_bending_and_stability(x, Sim, "GUI", "off", "DispInfo", "no", "InputAngleUnit", "Radian"), ...
    x_init, [], [], [], [], lb, ub, [], options);


fprintf('==== Optimization Finished ====\n');
fprintf('Best Cost (fval_opt): %.4f\n', f_val);
disp('Optimized x_opt = [theta1, theta2, ...]:');
disp(rad2deg(x_opt(1:6)));
toc;