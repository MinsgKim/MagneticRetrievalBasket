function obj = optimize_magnetic_robot_second()
    clc; clear; close all;

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

    % parameters of an external magnet
    mu0 = 4 * pi * 1e-7; % vacuum permeability
    external_magnet.Br = 1.22; % remanence [T]
    external_magnet.volume = (0.03)^3; % Volume (3cm x 3cm x 3cm cube)
    external_magnet.m = external_magnet.Br * external_magnet.volume / mu0; % magnetic moment
    external_magnet.position = [-0.00165; num_links * link_length_init + r_init]; % position of the magnet

    % optimized parameters storage
    num_iterations = 1; % iterations
    x_results = zeros(num_iterations, length(x0));
    cost_values = zeros(num_iterations, 1);

    % option setup
    options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-8, 'MaxFunctionEvaluations', 1e5, ...
        'FiniteDifferenceStepSize', 1e-6, 'OptimalityTolerance', 1e-8, 'Algorithm', 'interior-point');

    for i = 1:num_iterations
        tic;
        % input initial values
        x0 = [psi_init, theta_M_init, r_init, link_length_init];

        % set optim. problem
        problem = createOptimProblem('fmincon', 'x0', x0, ...
            'objective', @(x) objective_function(x, num_links, external_magnet), ...
            'lb', lb, 'ub', ub, ...
            'nonlcon', @(x) nonlcon_position_constraints(x, num_links, external_magnet), ...
            'options', options);

        % optimize!
        [x_result, fval] = fmincon(problem);

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
    [T_actual_opt, theta_opt] = simulate_robot_transform(num_links, M_opt, theta_M_opt, r_opt, link_length_opt, external_magnet);

    for k = 1:num_links
        obj(:, 3 * k - 2:3 * k) = T_actual_opt{k};
    end

    % robot angle
    disp(['theta (degrees): ', num2str(rad2deg(theta_opt))]);

    % robot visualization
    plot_robot(T_actual_opt);
end

function cost = objective_function(x, num_links, external_magnet)
    % cost function: maximize summation of middle links -> bending a lot
    % optim variables
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    r = x(2*num_links+1);
    link_length = x(end);
    cross_section_area = 0.0033 * 0.0005; % cross sectional area
    M = M * link_length * cross_section_area;

    % magnetic robot simulation
    [T_actual, ~] = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet);

    % calculate the sum of middle links
    x_positions = zeros(num_links-2, 1);
    for i = 2:num_links-1
        x_positions(i-1) = T_actual{i}(1, 3);
    end
    x_sum = sum(x_positions);

    % To maximize, mimimize the cost
    cost = -x_sum;
end

function [c, ceq] = nonlcon_position_constraints(x, num_links, external_magnet)
    % nonlinear constraints setup
    % optim variables
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    r = x(2*num_links+1);
    link_length = x(end);
    cross_section_area = 0.0033 * 0.0005; % cross sectional area
    M = M * link_length * cross_section_area;

    % magnetic robot simulation
    [T_actual, ~] = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet);

    % extract the position of each link
    positions = zeros(2, num_links);
    for i = 1:num_links
        positions(:, i) = T_actual{i}(1:2, 3);
    end

    % equation constraints: x position of 1st and last links are same
    ceq = positions(1, end) - positions(1, 1);

    % nonequation constraints:
    c = positions(1, 1) - positions(1, 2:end-1) + 1e-6; % minute tolerance 1e-6
    y_con = 0.005 - positions(2, end); % y coord. of last link > 5 mm (not negative)
    c = [c(:); y_con];
end

function [T_actual, theta_final] = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet)
    % parameter setup
    k_spring = 9.21e-4 * ones(1, num_links - 1); % spring coefficient of joint (PDMS or Ecoflex)
    damping_coefficient = 1e-6; % damping coefficient
    theta_init = 1e-3 * randn(1, num_links); % initial random angle of each link
    t_span = [0, 500]; % simulation time
    options_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);

    % simulation with ode15s
    [~, Y] = ode15s(@(t, y) robot_dynamics(t, y, k_spring, M, theta_M, r, link_length, external_magnet, damping_coefficient), ...
                   t_span, [theta_init, zeros(1, num_links)], options_ode);

    % extract the final link angle
    theta_final = Y(end, 1:num_links);

    % calculate transformation matrix
    T_actual = compute_transform_matrices(theta_final, link_length);
end


function dY = robot_dynamics(~, Y, k_spring, M, theta_M, r, L, external_magnet, damping_coefficient)
    % state vector Y = [theta1, theta2, ..., omega1, omega2, ...]
    num_links = length(M);
    theta = Y(1:num_links);
    omega = Y(num_links+1:end);

    % initiate acceleration vector
    alpha = zeros(num_links, 1);

    % inertia moment (unit mass, length L)
    I = (1/12) * L^2;

    % calculate the position of each link
    positions = compute_link_positions(theta, L);

    % calculate torque on each link
    for i = 1:num_links
        % vector from an external magnet to a link
        r_vec = positions(:, i) - [-0.00165; num_links * L + r];

        % calculate magnetic field on a link
        B_ext = calculate_B_vector(r_vec, external_magnet.m);

        % magnetic moment vector of a link
        theta_moment = theta(i) + theta_M(i);
        m_i = M(i) * [sin(theta_moment); cos(theta_moment); 0];

        % calculate magnetic torque (only consider z)
        tau_magnetic = m_i(2) * B_ext(1) - m_i(1) * B_ext(2);

        % spring torque
        tau_spring = 0;
        if i > 1
            tau_spring = tau_spring - k_spring(i - 1) * (theta(i) - theta(i - 1));
        end
        if i < num_links
            tau_spring = tau_spring + k_spring(i) * (theta(i + 1) - theta(i));
        end

        % damping torque
        tau_damping = -damping_coefficient * omega(i);

        % find angle acceleration
        alpha(i) = (tau_magnetic + tau_spring + tau_damping) / I;
    end

    dY = [omega; alpha];
end

function positions = compute_link_positions(theta, link_length)
    % calculate the end position of each link
    num_links = length(theta);
    x = zeros(1, num_links);
    y = zeros(1, num_links);
    current_theta = 0;
    x_current = 0;
    y_current = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        x_current = x_current + link_length * sin(current_theta);
        y_current = y_current + link_length * cos(current_theta);
        x(i) = x_current;
        y(i) = y_current;
    end

    positions = [x; y]; % end point of each link
end

function T = compute_transform_matrices(theta, link_length)
    % calculate 2d transformation matrix
    num_links = length(theta);
    T = cell(1, num_links);
    x = 0;
    y = 0;
    current_theta = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        dx = link_length * sin(current_theta);
        dy = link_length * cos(current_theta);

        % update position
        x = x + dx;
        y = y + dy;

        % 2d transformation matrix
        T{i} = [cos(current_theta), -sin(current_theta), x;
                sin(current_theta),  cos(current_theta), y;
                0,                  0,                  1];
    end
end

function plot_robot(T_actual)
    figure;
    hold on;
    num_links = length(T_actual);
    x_positions = zeros(1, num_links + 1);
    y_positions = zeros(1, num_links + 1);
    x_positions(1) = 0;
    y_positions(1) = 0;
    for i = 1:num_links
        x_positions(i + 1) = T_actual{i}(1, 3);
        y_positions(i + 1) = T_actual{i}(2, 3);
    end
    plot(x_positions, y_positions, '-o', 'LineWidth', 2);
    xlabel('X Position');
    ylabel('Y Position');
    title('Optimized Robot Configuration');
    grid on;
    axis equal;
    hold off;
end

function B = calculate_B_vector(r_vec, m_ext)
    % dipole model to calculate B-field
    mu0 = 4 * pi * 1e-7; % vacuum permeability

    r_norm = norm(r_vec);

    % magnetic moment of an external magnet (toward bottom)
    m_vec = [0; -m_ext; 0];

    % vector as 3D
    r_vec_3D = [r_vec; 0];

    % B-field
    B_full = (mu0 / (4 * pi)) * ((3 * r_vec_3D * (dot(m_vec, r_vec_3D)) / r_norm^5) - (m_vec / r_norm^3));

    % B-field in 2D (x, y components)
    B = B_full(1:2);
end
