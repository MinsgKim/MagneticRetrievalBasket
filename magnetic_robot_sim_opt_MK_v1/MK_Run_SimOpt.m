%% Preprocessing
clear; clc;
addpath("MatlabFunctions")
MK_Simscape_Preprocessing;

%%
Fitness = MK_Simscape_Fitness(Input.x_init_1, SimParams, Fixed, "DispInfo", "yes", "InputAngleUnit", "Radian", "GUI", "on");
% Plot_Config(x, )

%% Optimizing Simulation

lb = [repmat(-pi, 1, Fixed.num_links), repmat(5e3, 1, Fixed.num_links) * Fixed.PlateVolume];
ub = [repmat(pi, 1, Fixed.num_links), repmat(5e4, 1, Fixed.num_links) * Fixed.PlateVolume];

options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-14, 'FunctionTolerance', 1e-14, ...
    'MaxFunctionEvaluations', 1e5, 'OptimalityTolerance', 1e-7, ...
    'Algorithm', 'sqp');

problem = createOptimProblem('fmincon',...
    'x0', Input.x_init_random, ...
    'lb', lb, 'ub', ub, ...
    'objective', @(x) Simscape_Fitness(x, SimParams, Fixed),...
    'options', options);

gs = GlobalSearch;

tic;

num_iterations = 1;

for i=1:num_iterations
    [x_opt, fval_opt] = run(gs, problem);
    x_results(i,:) = x_opt;
    cost_values(i) = fval_opt;
    fprintf("Iteration %d completed.\n", i)
end

[~, best_cost] = min(cost_values);
x_opt = x_results(best_cost, :);

% [x_opt, f_val] = fmincon(@(x) Simscape_Fitness(x, SimParams, Fixed), ...
%     Input.x_init_random, [], [], [], [], lb, ub, [], options);

fprintf('==== Optimization Finished ====\n');
fprintf('Best Cost (fval_opt): %.4f\n', f_val);
disp('Optimized x_opt = [thetaM_1..thetaM_n]:');
disp(rad2deg(x_opt(1:Fixed.num_links)));
disp('Optimized x_opt = [psi_1..psi_n]:');
disp(x_opt(Fixed.num_links+1:end)/Fixed.PlateVolume);
toc;

%% Calibration of the permanent magnet
%-----------objective: find a dipole moment of a single magnet------------%
clc; clear;

filename = 'cubic_magnetic_field_data2.xlsx';

sheetname = 'sheet1';
% range = 'A1:C100';

dataTable = readtable(filename, 'Sheet', sheetname); % data step 1 mm??

doubleData = table2array(dataTable);

normData = doubleData(:, end) * 1e-3;

for i = 1:length(normData)

    B(:,1,i) = [normData(i); 0; 0];
    % B(:,1,i) = [-doubleData(i,1); -doubleData(i,2); doubleData(i,3)] * 1e-3;
    r(:,1,i) = [(5 + 40.5 + 1*(i-1)); 0; 0] * 1e-3;
    % r(:,1,i) = [(5 + 40.5 + 1*(i-1)); -1; -1] * 1e-3;
    r_hat(:,1,i) = r(:,1,i)/pagenorm(r(:,1,i));

    mu = 4*pi*1e-7;

    A(:,:,i) = mu/(4*pi*pagenorm(r(:,1,i))^3)*(3*(r_hat(:,1,i)*r_hat(:,1,i)')-eye(3,3));

    M(i) = (inv(A(:,:,i))*B(:,1,i))'*[1; 0; 0];

end

t = 1:length(M);

plot(t, M)

values = [min(M), max(M(2:end)), mean(M(2:end))];
disp(values)

%%

% B_small_1 = [-0.072391; 0.1129425; 0.02058] * 1e-3;
% B_small_2 = [-0.0318769; -0.104541; -0.0023366] * 1e-3;

B = [-0.04; -0.01; -0.01]*1e-3;

r = [10.5+0.3; 10; 0] * 1e-3;
r_hat = r/norm(r);

mu = 4*pi*1e-7;

A = mu/(4*pi*norm(r)^3)*(3*(r_hat*r_hat')-eye(3,3));

M_1 = inv(A) * B_small_1;
M_2 = inv(A) * B_small_2;

M_norm_1 = norm(M_1);
M_norm_2 = norm(M_2);

M = inv(A) * B;

M_norm = norm(M);
disp(M_norm)