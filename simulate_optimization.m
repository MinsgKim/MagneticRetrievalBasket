clc; clear; close all;

EM = External_Magnet;
mrs = magnetic_robot_simulation;
RK = Robot_Kinematics;


% optimizing parameters
num_links = 7; % the number of links
psi_init = 4e04 * ones(1, num_links); % initial magnetization profile [A/m]
rng(0); % fix random generator
theta_M_init = rand(1, num_links) * 2 * pi - pi; % magnetization direction initial values (0)
r_init = 0.045; % initial distance from an external magnet to the robot end [m]
link_length_init = 2e-03; % link length

obj = zeros(3, 3 * num_links);

% initial guess
x0 = [psi_init, theta_M_init, r_init, link_length_init];

% optimizing boundaries
lb = [repmat(2e04, 1, num_links), repmat(-pi, 1, num_links), 0.02, 0.001];
ub = [repmat(6e04, 1, num_links), repmat(pi, 1, num_links), 0.05, 0.003];


% optimized parameters storage
num_iterations = 5; % iterations
x_results = zeros(num_iterations, length(x0));
cost_values = zeros(num_iterations, 1);

% option setup
options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6, 'MaxFunctionEvaluations', 1e5, ...
    'FiniteDifferenceStepSize', 1e-6, 'OptimalityTolerance', 1e-6, 'Algorithm', 'interior-point',"EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg");

for i = 1:num_iterations
    tic;
    % input initial values
    x0 = [psi_init, theta_M_init, r_init, link_length_init];

    % set optim. problem
    problem = createOptimProblem('fmincon', 'x0', x0, ...
        'objective', @(x) mrs.objective_function(x, num_links), ...
        'lb', lb, 'ub', ub, ...
        'nonlcon', @(x) mrs.nonlcon_position_constraints(x, num_links), ...
        'options', options);

    % optimize!
    gs = GlobalSearch;
    [x_result, fval] = run(gs, problem);

    x_results(i, :) = x_result;
    cost_values(i) = fval;
    disp(['Iteration ', num2str(i), ' completed.']);
    toc;
end



% choose values with a minimal cost
[~, best_idx] = min(cost_values);
x_opt = x_results(best_idx, :);

% display param
disp('optimized parameters:');
disp(['M: ', num2str(x_opt(1:num_links))]);
disp(['theta_M (degrees): ', num2str(rad2deg(x_opt(num_links+1:2*num_links)))]);
disp(['r: ', num2str(x_opt(end-1))]);
disp(['link_length: ', num2str(x_opt(end))]);

% robot simulation & visualization
M_opt = x_opt(1:num_links);
theta_M_opt = x_opt(num_links+1:2*num_links);
r_opt = x_opt(end-1);
link_length_opt = x_opt(end);
cross_section_area = 0.0033 * 0.0005; % cross sectional area (3.3 mm x 0.5 mm)
M_opt = M_opt * link_length_opt * cross_section_area;

% robot simulation
[T_actual_opt, theta_opt] = RK.simulate_robot_transform(num_links, M_opt, theta_M_opt, r_opt, link_length_opt, EM);

for k = 1:num_links
    obj(:, 3 * k - 2:3 * k) = T_actual_opt{k};
end

% robot angle
disp(['theta (degrees): ', num2str(rad2deg(theta_opt))]);

% robot visualization
RK.plot_robot(T_actual_opt);