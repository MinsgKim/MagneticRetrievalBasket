% Let's set a modified optimization problem.
% This time we are going to lay down a robot so that we can start an angle
% from x axis.
% An external magnet will be placed on the right side of the coordinates.
% Just a different POV. It is same when we rotate the screen 90 degrees.

% initial parameter setup
clear; clc; close all

num_links = 7;  % # of links
link_length = 2e-3; % 2 mm length of each link

psi_init = 5000 * (7 * rand(1,num_links) + 1);          % magnetization profile [A/m]
theta_M_init = pi/4 * (2* rand(1, num_links) - 1);      % magnetization angle [rad]
theta_link_init = pi/16 * (2 * rand(1, num_links) - 1);  % link angle [rad]

r_ext = 30e-3;  % from center of an external magnet to the first joint in x direction
cross_section_area = 0.0033*0.0005;     % cross section area of each part
k_spring = 1e-05 * ones(1,7);            % spring constants of all joints

% build class objects
cf = cost_function();
em2 = External_Magnet2();
RS = RobotState();

% initial guess
x_init = [psi_init, theta_M_init, theta_link_init];

% optimizing boundaries
lb = [repmat(5e03, 1, num_links), repmat(-pi, 1, num_links), repmat(-pi/4, 1, num_links)];
ub = [repmat(6e04, 1, num_links), repmat(pi, 1, num_links), repmat(pi/4, 1, num_links)];

% optimizing linear inequality
% A = [zeros(num_links, 2 * num_links), eye(num_links)
%     zeros(num_links, 2 * num_links), -eye(num_links)];
% b = [zeros(num_links * 2) ]';

% optimized parameters storage
num_iterations = 1; % iterations
x_results = zeros(num_iterations, length(x_init));
cost_values = zeros(num_iterations, 1);

% option setup
options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6, 'MaxFunctionEvaluations', 1e5, ...
    'OptimalityTolerance', 1e-6, 'Algorithm', 'interior-point',"EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg");
% createOptimProblem + GlobalSearch
problem = createOptimProblem('fmincon',...
    'x0', x_init, ...
    'lb', lb, 'ub', ub, ...
    'objective', @(x) cf.moment_equilibrium(x, num_links, link_length, cross_section_area, r_ext, k_spring, em2),...
    'nonlcon', @(x) cf.nonlcon(x, num_links, link_length, cross_section_area, r_ext, k_spring, em2),...
    'options', options);

gs = GlobalSearch;
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
disp('Optimized x_opt = [theta_1..theta_n]:');
disp(rad2deg(x_opt(2*num_links+1:end)));

cost_examination = cf.moment_equilibrium(x_opt, num_links, link_length, cross_section_area, r_ext, k_spring, em2);
T = cf.get_tau(x_opt, num_links, link_length, cross_section_area, r_ext, k_spring, em2);

link_center_opt = RS.set_link_center(num_links, link_length, x_opt(2*num_links:end));

from_mag_to_link = zeros(2, num_links);
B_field = zeros(3, num_links);

for i = 1:num_links

    from_mag_to_link(:, i) = link_center_opt(:,i) - [r_ext; -0.00165];
    B_field(:, i) = em2.Cal_B_Field(from_mag_to_link(:, i));

end


RS.draw_plot(num_links, link_length, x_opt(2*num_links+1:end), x_opt(num_links+1:2*num_links))


%% function examination
clear; clc; close all

num_links = 7;  % # of links
link_length = 2e-3; % 2 mm length of each link

psi_init = [40000 20000 10000 10000 10000 10000 10000];          % magnetization profile [A/m]
theta_M_init = [-pi/2 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4];      % magnetization angle [rad]
theta_link_init = [pi/4 -pi/16 -pi/16 -pi/16 -pi/16 -pi/16 -pi/16];  % link angle [rad]

r_ext = 35e-3;  % from center of an external magnet to the first joint in x direction
cross_section_area = 0.0033*0.0005;     % cross section area of each part
k_spring = 3e-05 * ones(1,7);            % spring constants of all joints

% build class objects
cf = cost_function();
em2 = External_Magnet2();
RS = RobotState();

% test
x_test = [psi_init, theta_M_init, theta_link_init];

cost = cf.moment_equilibrium(x_test, num_links, link_length, cross_section_area, r_ext, k_spring, em2);

link_center_test = RS.set_link_center(num_links, link_length, x_test(2*num_links:end));

for i = 1:num_links

    from_mag_to_link = link_center_test(:,i) - [r_ext; -0.00165];
    B_field = em2.Cal_B_Field(from_mag_to_link);
    fprintf('distance: %.4e\n B field: %.4e\n', from_mag_to_link, norm(B_field))
    pause(0.1)

end


RS.draw_plot(num_links, link_length, x_test(2*num_links+1:end), x_test(num_links+1:2*num_links))


%% really really simple test

clear;
clc;

em1 = External_Magnet;
em2 = External_Magnet2;

r1 = [6.12323399573677e-20	0.000217392223869186	0.000791492065212027	0.00153735007829041	0.00224904192297457	0.00265898082881719	0.00231199527224396
-0.00100000000000000	-0.00297608433088602	-0.00488638473722582	-0.00674177504679641	-0.00860950462193969	-0.0105522334842722	-0.0124491352652242];

r2 = [r1(2,:); r1(1,:)];

r3 = repmat([35e-3; -0.00165], 1, 7);
r4 = repmat([-0.00165; 35e-03], 1, 7);

r_x = r1-r3;
r_y = r2-r4;

b_field1 = [];
b_field2 = [];
for i = 1:7
    b_field1(i) = norm(em1.Cal_B(r_y(:,i)));
    b_field2(i) = norm(em2.Cal_B_Field(r_x(:,i)));

end


% distance: -3.5000e-02
%  B field: 6.5000e-04
% distance: 1.6360e-02
%  B field: distance: -3.4783e-02
%  B field: -1.3261e-03
% distance: 1.6634e-02
%  B field: distance: -3.4209e-02
%  B field: -3.2364e-03
% distance: 1.7243e-02
%  B field: distance: -3.3463e-02
%  B field: -5.0918e-03
% distance: 1.7945e-02
%  B field: distance: -3.2751e-02
%  B field: -6.9595e-03
% distance: 1.8393e-02
%  B field: distance: -3.2341e-02
%  B field: -8.9022e-03
% distance: 1.8098e-02
%  B field: distance: -3.2688e-02
%  B field: -1.0799e-02
% distance: 1.6556e-02
%% just...just test...

clear; clc;

em2 = External_Magnet2();
M = 45000 * 0.0033 * 0.0005 * 2e-03;
theta_M = pi/2;
link_length = 2e-03;
r_ext = [35e-03; -0.00165];
theta_angle = linspace(0, pi/3, 10);

k = 2.51e-6;

B_Field = zeros(3, length(theta_angle));
M_vec = zeros(3, length(theta_angle));
t_mag = zeros(length(theta_angle),1);
t_k = zeros(length(theta_angle),1);

diff = zeros(length(theta_angle),1);

for i = 1:length(theta_angle)

    link_center = link_length/2*[cos(theta_angle(i)); sin(theta_angle(i))];
    r_vec = link_center - r_ext;
    B_Field(:,i) = em2.Cal_B_Field(r_vec);
    M_vec(:, i) = M * [cos(theta_M + theta_angle(i)); sin(theta_M + theta_angle(i)); 0];
    t_mag_vec = cross(M_vec(:,i), B_Field(:,i));
    t_mag(i) = t_mag_vec(3);
    t_k(i) = k * theta_angle(i);
    diff(i) = t_mag(i) - t_k(i);

end

t_mag = t_mag * 1e+7;
t_k = t_k * 1e+7;
diff = diff *1e+7;

disp(rad2deg(theta_angle))

% disp(num2str(t_mag))
% disp(num2str(t_k))
% disp(num2str(diff))


%% operator version-- 2








