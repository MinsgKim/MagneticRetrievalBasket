%% optimize for calibration
clear;clc;
Sim = Call_MagneticBasket_Params;
open("Calibration_250328.slx")

Sim.Angle = deg2rad([25 30 45 -45 -35 -25]);
% Sim.Transition_upper = 0.1745; % 10 degree
Sim.Transition_upper = 0.6; %

% data
measured_data(:,:,1) = [738 136; 775 130; 811 126; 852 120; 888 107;...
    923 80; 943 47];
measured_data(:,:,2) = [734 138; 774 131; 811 127; 852 122; 882 92;...
    897 53; 901 11];
measured_data(:,:,3) = [734 138; 771 133; 810 129; 853 123; 879 87;...
    888 48; 881 8];
measured_data(:,:,4) = [729 139; 774 137; 813 142; 853 153; 892 159;...
    932 166; 971 166];
measured_data(:,:,5) = [730 140; 775 139; 814 148; 849 160; 891 172;...
    929 185; 967 194];
measured_data(:,:,6) = [734 139; 777 148; 810 169; 839 202; 862 232;...
    886 268; 908 301];

measured_angle = zeros(1,7);
measured_q = zeros(1,6,6);

for j=1:6
    for i=1:6

        x(i) = measured_data(i,1,j)-measured_data(i+1,1,j);
        y(i) = measured_data(i,2,j)-measured_data(i+1,2,j);
        measured_angle(i+1) = rad2deg(atan(abs(y(i)./x(i))));
        measured_q(1,i,j) = measured_angle(i+1) - measured_angle(i);
    end
end

measured_q(:,:,4:6) = -measured_q(:,:,4:6);
measured_q = deg2rad(measured_q);

% disp(measured_q)
% disp(deg2rad(measured_q))

x_init = [[30973 30962 39510 25856 25211 25256]*Sim.PlateVolume, [2.5147e-6 2.2831e-6]* 1e2] * 1e4;

lb = [ones(1,6)*25000*Sim.PlateVolume, [0.1e-6 0.1e-6]* 1e2]*1e4;
ub = [ones(1,6)*40000*Sim.PlateVolume, [10e-6 10e-6] * 1e2]*1e4;


% x_init = [[8 8 8 8 8 11]*1e-4, [7e-5 8e-5]];
% 
% lb = [[7 7 7 7 7 7]*1e-4, [0.1e-5 0.1e-5]];
% ub = [[13 13 13 13 13 13]*1e-4, [10e-5 10e-5]];


options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-3, 'FunctionTolerance', 1e-3, ...
    'MaxFunctionEvaluations', 1e3, 'OptimalityTolerance', 1e-3, ...
    'Algorithm', 'interior-point', "EnableFeasibilityMode", true, ...
    "SubproblemAlgorithm","cg");

tic;
[x_opt, f_val] = fmincon(@(x) calibrate_Kt(x, measured_q, Sim, "GUI", "off", "DispInfo", "no"), ...
    x_init, [], [], [], [], lb, ub, [], options);


fprintf('==== Optimization Finished ====\n');
fprintf('Best Cost (fval_opt): %.4f\n', f_val);
disp('Optimized x_opt = [m1, m2, ...]:');
disp(x_opt(1:6));
disp('Optimized stiffness = Kt:');
disp(x_opt(end));
toc;