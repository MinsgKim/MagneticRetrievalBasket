%% operator version-- 2
% you should run section

clear; clc; close all

num_links = 7;  % # of links
link_length = 2e-3; % 2 mm length of each link

psi_init = 5000 * (7 * rand(1,num_links) + 1);          % magnetization profile [A/m]
theta_M_init = [pi/3, pi/4 * (2* rand(1, num_links-1) - 1)];      % magnetization angle [rad]

r_ext = 45e-3;  % from center of an external magnet to the first joint in x direction
cross_section_area = 0.003*0.0005;     % cross section area of each part
k_spring = 2.5e-06 * ones(1,7);            % spring constants of all joints

% build class objects
cf = cost_function();
em2 = External_Magnet2();
RS = RobotState();

% initial guess
x_init = [psi_init, theta_M_init];

% optimizing boundaries
lb = [[40000, repmat(5e03, 1, num_links-1)], [pi/3, repmat(-pi, 1, num_links-1)]];
ub = [[50000, repmat(5e04, 1, num_links-1)], [pi/2, repmat(pi, 1, num_links-2), 0]];

% optimized parameters storage
num_iterations = 2; % iterations
x_results = zeros(num_iterations, length(x_init));
cost_values = zeros(num_iterations, 1);

% option setup
options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-5, 'FunctionTolerance', 1e-5, ...
    'ConstraintTolerance', 5e-5, 'MaxFunctionEvaluations', 1e5, ...
    'OptimalityTolerance', 1e-5, 'Algorithm', 'interior-point', "EnableFeasibilityMode", true,...
    "SubproblemAlgorithm","cg");

% createOptimProblem + GlobalSearch
problem = createOptimProblem('fmincon',...
    'x0', x_init, ...
    'lb', lb, 'ub', ub, ...
    'objective', @(x) cf.Max_Bending(x, num_links, link_length, cross_section_area, r_ext, k_spring, em2),...
    'nonlcon', @(x) cf.nonlcon_ver2(x, num_links, link_length, cross_section_area, r_ext, k_spring, em2),...
    'options', options);

gs = GlobalSearch;
tic;
for i=1:num_iterations
    [x_opt, fval_opt] = run(gs, problem);
    x_results(i,:) = x_opt;
    cost_values(i) = fval_opt;
    fprintf("Iteration %d completed.\n", i)
end

[~, best_cost] = min(cost_values);
x_opt = x_results(best_cost, :);

fprintf('==== Optimization Finished ====\n');
fprintf('Best Cost (fval_opt): %.4f\n', min(cost_values));
disp('Optimized x_opt = [psi_1..psi_n]:');
disp(x_opt(1:num_links));
disp('Optimized x_opt = [thetaM_1..thetaM_n]:');
disp(rad2deg(x_opt(num_links+1:2*num_links)));

cost_examination = cf.Max_Bending(x_opt, num_links, link_length, cross_section_area, r_ext, k_spring, em2);
M_opt = x_opt(1:num_links) * cross_section_area * link_length;
theta_opt = RS.Get_Link_Angle(num_links, link_length, M_opt, x_opt(num_links+1:end), r_ext, k_spring, em2);
disp(rad2deg(theta_opt))

[T_m, T_s, T_sum] = RS.Get_Tau(num_links, link_length, M_opt, x_opt(num_links+1:end), r_ext, k_spring, em2);

% from_mag_to_link = zeros(2, num_links);
% B_field = zeros(3, num_links);


RS.draw_plot(num_links, link_length, theta_opt, x_opt(num_links+1:end))
toc;

%% forward test version --2


clear; clc; close all

num_links = 7;  % # of links
link_length = 2e-3; % 2 mm length of each link

% psi_init = [44000 31000 16000 48000 24000 50000 40000];         % magnetization profile [A/m]
psi_init = [44000 31000 16000 48000 24000 50000 40000];
theta_M_init = deg2rad([75 62 -40 -100 -95 -91 -80]);      % magnetization angle [rad]
% theta_M_init = [pi/2 pi/2 pi/2 -pi/2 -pi/2 -pi/2 -pi/2];

r_ext = 45e-3;  % from center of an external magnet to the first joint in x direction
cross_section_area = 0.003*0.0005;     % cross section area of each part
k_spring = 2.5e-06 * ones(1,num_links);            % spring constants of all joints

% build class objects
cf = cost_function();
em2 = External_Magnet2();
RS = RobotState();

x_init = [psi_init, theta_M_init];

M_opt = x_init(1:num_links) * cross_section_area * link_length;
theta_opt = RS.Get_Link_Angle(num_links, link_length, M_opt, x_init(num_links+1:end), r_ext, k_spring, em2);

[T_m, T_s, T_sum, isConverged] = RS.Get_Tau(num_links, link_length, M_opt, x_init(num_links+1:end), r_ext, k_spring, em2);

RS.draw_plot(num_links, link_length, theta_opt, x_init(num_links+1:end))

%%

em2 = External_Magnet2;

r = [0.044; -0.001587];

B = em2.Cal_B_Field(r);

