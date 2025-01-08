clc; clear; close all;

EM = External_Magnet;
mrs = magnetic_robot_simulation;
RK = Robot_Kinematics;


% optimizing parameters
num_links = 7; % the number of links
psi_init = 5000 * (7 * rand(1, num_links) + 1.0); % initial magnetization profile [A/m]
% rng(0); % fix random generator
theta_M_init = rand(1, num_links) * 2 * pi - pi; % magnetization direction initial values (0)
r_init = 0.05; % initial distance from an external magnet to the robot end [m]
link_length_init = 2e-03; % link length

obj = zeros(3, 3 * num_links);

% initial guess
x0 = [psi_init, theta_M_init, r_init, link_length_init];

% optimizing boundaries
lb = [repmat(5e03, 1, num_links), repmat(-pi, 1, num_links), 0.04, 0.001];
ub = [repmat(4e04, 1, num_links), repmat(pi, 1, num_links), 0.06, 0.003];


% optimized parameters storage
num_iterations = 1; % iterations
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
RK.plot_robot(T_actual_opt, theta_opt, theta_M_opt, link_length_opt, r_opt);

%% simple test

clc; clear; close all;

EM = External_Magnet;
mrs = magnetic_robot_simulation;
RK = Robot_Kinematics;

% optimizing parameters
num_links = 12; % the number of links
% psi_init = 1e03 * (rand(1, num_links) + 0.5); % initial magnetization profile [A/m]
% psi_init = 1e04 * ones(1, num_links);   % 30,000 A/m -> 1:1 ratio
psi_init = [60000 40000 20000 20000 20000 10000 10000 10000 10000 10000 10000 10000];
rng(0); % fix random generator
% theta_M_init = -[-pi -pi -pi -pi pi pi pi]/2; % magnetization direction initial values (0)
theta_M_init = -[-pi/2 -pi/2 -pi/2 -pi/2 -pi/2 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/2]; % magnetization direction initial values (0)
r_init = 0.03; % initial distance from an external magnet to the robot end [m]
link_length_init = 1e-03; % link length
cross_section_area = 0.0033 * 0.0005;
M_init = psi_init * link_length_init * cross_section_area;

x0 = [psi_init, theta_M_init, r_init, link_length_init];

cost = mrs.objective_function(x0, num_links);

[T_test, theta_test] = RK.simulate_robot_transform(num_links, M_init, theta_M_init, r_init, link_length_init, EM);

RK.num_links = num_links;
pos = RK.compute_link_positions(theta_test, link_length_init);

for i=1:num_links

    if i > 1
        x_ct = (pos(1,i)+pos(1,i-1))/2;
        y_ct = (pos(2,i)+pos(2,i-1))/2;

    else
        x_ct = pos(1,i)/2;
        y_ct = pos(2,i)/2;

    end

    r_vec = [x_ct; y_ct]-[0.00165; r_init];
    fprintf('%d번째:\n', i)
    disp(['distance: ', num2str(norm(r_vec))])
    disp(['B-field: ', num2str(norm(EM.Cal_B(r_vec)))])

end

RK.plot_robot(T_test, theta_test, theta_M_init, link_length_init, r_init);

%% static test

clc; clear; close all;

%%% 전역 변수 선언 (최적화 중 기록용) %%%
global best_fval best_theta
best_fval = inf;
best_theta = [];

% 로봇 링크 개수
num_links = 7;
link_length = 2e-3;
cross_section_area = 0.0033 * 0.0005;
k_spring = 3e-5 * ones(1, num_links-1);
r_ext = 0.04;
EM = External_Magnet();
mrs = magnetic_robot_simulation();
RK = Robot_Kinematics();

% ---- 최적화 변수 정의 (psi, theta_M) ----
lb_psi = 5e3;  ub_psi = 6e4;
lb_thetaM = -pi;  ub_thetaM = pi;
lb = [lb_psi*ones(1,num_links), lb_thetaM*ones(1,num_links)];
ub = [ub_psi*ones(1,num_links), ub_thetaM*ones(1,num_links)];

% 초기값
psi_init = 5e3 * (7 * rand(1,num_links) + 1);
thetaM_init = (rand(1,num_links)*2*pi - pi);
x0 = [psi_init, thetaM_init];

% optimized parameters storage
num_iterations = 1;
x_results = zeros(num_iterations, length(x0));
cost_values = zeros(num_iterations, 1);

% fmincon 옵션
options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10, 'MaxFunctionEvaluations', 1e5, ...
    'FiniteDifferenceStepSize', 1e-6, 'OptimalityTolerance', 1e-6, 'Algorithm', 'interior-point',"EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg");

% createOptimProblem + GlobalSearch
problem = createOptimProblem('fmincon',...
    'x0', x0, ...
    'lb', lb, 'ub', ub, ...
    'objective', @(x) mrs.objective_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
    'nonlcon', @(x) mrs.constraint_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
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
disp('Optimized x_opt = [psi_1..psi_n, thetaM_1..thetaM_n]:');
disp(x_opt);

%%% 여기서 더 이상 solve_static_equilibrium을 새로 부르지 않고,
%%% "best_theta"를 최종 로봇 형상으로 활용
%%% --------------------------------------------
global best_theta
if isempty(best_theta)
    warning('best_theta is empty! Possibly no update occurred?');
    return;
end

theta_final = best_theta;
positions_opt = RK.compute_link_positions2(theta_final, link_length);

% 검산: 중간 링크 x좌표 합
if num_links>2
    x_mid = positions_opt(1,2:end-1);
else
    x_mid = 0;
end
cost_examination = -sum(x_mid);
fprintf('Re-check cost with best_theta = %.4f\n', cost_examination);

%%% 최종 플롯
% 자기장 계산
figure; hold on; axis equal;

for i = 1:num_links
    r_vec = positions_opt(:,i) - [0; r_ext];
    B_local = EM.Cal_B(r_vec);
    fprintf('Link %d:  r=%.4e, |B|=%.4e\n', i, norm(r_vec), norm(B_local));
end

% 로봇 형상
for i = 1:num_links-1
    plot(positions_opt(1,i:i+1), positions_opt(2,i:i+1), 'bo-','LineWidth',2 );
end
plot([0 positions_opt(1,1)], [0 positions_opt(2,1)], 'bo-','LineWidth',2)
plot([0, positions_opt(1,:)], [0, positions_opt(2,:)], 'ro','MarkerSize',8,'LineWidth',2 );

% 자화 방향 화살표
psi_opt    = x_opt(1:num_links);
thetaM_opt = x_opt(num_links+1 : 2*num_links);
current_th = 0;
for i = 1:num_links
    current_th = current_th + theta_final(i);
    th_m = theta_final(i) + thetaM_opt(i);
    if i == 1
        x_center = positions_opt(1,i)/2;
        y_center = positions_opt(2,i)/2;
    else
        x_center = (positions_opt(1,i)+positions_opt(1,i-1))/2;
        y_center = (positions_opt(2,i)+positions_opt(2,i-1))/2;
    end
    quiver( x_center, y_center, 0.5*link_length*sin(th_m), 0.5*link_length*cos(th_m), ...
        'Color',[1,0,0], 'LineWidth',1.5, 'MaxHeadSize',2 );
end

title('Final Robot Configuration (Best Cost)');
xlabel('X (m)'); ylabel('Y (m)');
grid on;
hold off;

%% forward test (static)
clear; clc; close all;


num_links = 7;                              % 링크 개수
link_length = 2e-3;                         % 각 링크 길이
cross_section_area = 0.0033 * 0.0005;
k_spring = 3e-5 * ones(1, num_links-1);
r_ext = 0.04;                               % 외부 자석 거리
EM = External_Magnet();
mrs = magnetic_robot_simulation();
RK = Robot_Kinematics();

% 초기값
% psi_init = 5e3 * (7 * rand(1,num_links) + 1);
psi_init = [40000 10000 10000 5000 5000 5000 2000];
thetaM_init = [-pi/4 -pi/4 pi/4 pi/4 pi/4 pi/4 pi/4];
x0 = [psi_init, thetaM_init];

global theta_test2 

cost = mrs.objective_static(x0, num_links, link_length, cross_section_area, r_ext, k_spring, EM);

positions = RK.compute_link_positions2(theta_test2, link_length);
% 검산: 중간 링크 x좌표 합
if num_links>2
    x_mid = positions(1,2:end-1);
else
    x_mid = 0;
end
cost_examination = -sum(x_mid);

% 자기장 계산
figure; hold on; axis equal;

for i = 1:num_links
    r_vec = positions(:,i) - [0; r_ext];
    B_local = EM.Cal_B(r_vec);
    fprintf('Link %d:  r=%.4e, |B|=%.4e\n', i, norm(r_vec), norm(B_local));
end


% 로봇 형상
for i = 1:num_links-1
    plot(positions(1,i:i+1), positions(2,i:i+1), 'bo-','LineWidth',2 );
end
plot([0 positions(1,1)], [0 positions(2,1)], 'bo-','LineWidth',2)
plot([0, positions(1,:)], [0, positions(2,:)], 'ro','MarkerSize',8,'LineWidth',2 );

% 자화 방향 화살표
current_th = 0;
for i = 1:num_links
    current_th = current_th + theta_test2(i);
    th_m = theta_test2(i) + thetaM_init(i);
    if i == 1
        x_center = positions(1,i)/2;
        y_center = positions(2,i)/2;
    else
        x_center = (positions(1,i)+positions(1,i-1))/2;
        y_center = (positions(2,i)+positions(2,i-1))/2;
    end
    quiver( x_center, y_center, 0.5*link_length*sin(th_m), 0.5*link_length*cos(th_m), ...
        'Color',[1,0,0], 'LineWidth',1.5, 'MaxHeadSize',2 );
end

title('Final Robot Configuration (Best Cost)');
xlabel('X (m)'); ylabel('Y (m)');
grid on;
hold off;

