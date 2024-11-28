function optimize_magnetic_robot()
    % Optimization parameters
    num_links = 7; % Fixed number of links
    psi_init = 3e04*ones(1, num_links); % magnetization profiles (A/m), Initialize with all ones
    theta_M_init = rand(1, num_links); % magnetization direction (rad), Initialize with all zeros
    B_init = [0.01, pi/4]; % Magnitude and angle
    link_length_init = 2e-03; % 2 mm initial link length
    disp(rad2deg(theta_M_init))

    % Initial guess
    x0 = [psi_init, theta_M_init, B_init, link_length_init];

    % Optimization bounds
    lb = [repmat(2e04, 1, num_links), repmat(-pi, 1, num_links), [0, pi/4], 0.001];
    ub = [repmat(6e04, 1, num_links), repmat(pi, 1, num_links), [0.1, pi/2], 0.01];

    % Fixed position constraints
    x_fixed_last = 1e-10; % Fixed x position of the last link (example)

    % Optimization options
    options = optimoptions('fmincon','Display','iter','StepTolerance',1e-15,...
    'ConstraintTolerance',1e-9,'MaxFunctionEvaluations',8e3,...
    'FiniteDifferenceStepSize',1.5e-6,...
    'OptimalityTolerance',1e-8,'Algorithm','sqp');

    % Perform optimization
    x_opt = fmincon(@(x) objective_function_mid_position(x, num_links), ...
                    x0, [], [], [], [], lb, ub, ...
                    @(x) nonlcon_position_constraints(x, num_links, x_fixed_last), ...
                    options);

    % Display optimized parameters
    disp('Optimized Parameters:');
    disp(['M: ', num2str(x_opt(1:num_links))]);
    disp(['theta_M: ', num2str(rad2deg(x_opt(num_links+1:2*num_links)))]);
    disp(['B: ', num2str(x_opt(end-2:end-1))]);
    disp(['link_length: ', num2str(x_opt(end))]);
end

function cost = objective_function_mid_position(x, num_links)
    % Unpack optimization variables
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    B = x(end-2:end-1);
    link_length = x(end);
    M = M*x(end)*0.0033*0.0005;

    % External magnetic field in vector form
    B_vector = B(1) * [cos(B(2)); sin(B(2))];

    % Simulate the magnetic robot
    T_actual = simulate_robot_transform(num_links, M, theta_M, B_vector, link_length);

    % Extract position of the middle transformation matrix (3rd link)
    middle_index = ceil(num_links / 2);
    position = T_actual{middle_index}(1:2, 3); % Extract 2D position (x, y)
    position_norm = norm(position); % Calculate norm of the position

    % Maximize position norm by minimizing its negative
    cost = -position_norm;
end

function [c, ceq] = nonlcon_position_constraints(x, num_links, x_fixed_last)
    % Unpack optimization variables
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    B = x(end-2:end-1);
    link_length = x(end);
    M = M*x(end)*0.0033*0.0005;

    % External magnetic field in vector form
    B_vector = B(1) * [cos(B(2)); sin(B(2))];

    % Simulate the magnetic robot
    T_actual = simulate_robot_transform(num_links, M, theta_M, B_vector, link_length);

    % Constraint 1: First link's position is fixed (already at [0, 0])
    % No inequality constraints here as it's inherently satisfied

    % Constraint 2: Last link's x position is fixed
    last_position = T_actual{end}(1:2, 3); % Extract position of the last link
    ceq = last_position(1) - x_fixed_last; % x position must equal x_fixed_last

    % No inequality constraints
    c = [];
end

function T_actual = simulate_robot_transform(num_links, M, theta_M, B, link_length)
    % Parameters
    k_spring = 0.00092125 * ones(1, num_links-1);
    theta_init = zeros(1, num_links); % Initial link angles
    t_span = [0, 10];
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);

    % Simulate using ode45
    [~, Y] = ode45(@(t, y) robot_dynamics(t, y, k_spring, M, theta_M, B, link_length), ...
                   t_span, [theta_init, zeros(1, num_links)], options);

    % Compute transformation matrices
    T_actual = compute_transform_matrices(Y(end, 1:num_links), link_length);
end

function dY = robot_dynamics(t, Y, k_spring, M, theta_M, B, L)
    % State vector Y = [theta1, theta2, theta3, omega1, omega2, omega3]
    num_links = length(M);
    theta = Y(1:num_links);
    omega = Y(num_links+1:end);
    
    % Initialize acceleration vector
    alpha = zeros(num_links, 1);
    
    % Calculate torques for each link
    for i = 1:num_links
        % Magnetic torque
        tau_magnetic = M(i) * norm(B) * sin(theta(i) + theta_M(i) - atan2(B(2), B(1)));
        
        % Spring torques
        tau_spring = 0;
        if i > 1
            tau_spring = tau_spring - k_spring(i-1) * (theta(i) - theta(i-1));
        end
        if i < num_links
            tau_spring = tau_spring + k_spring(i) * (theta(i+1) - theta(i));
        end
        
        % Simple damping
        tau_damping = -0.01 * omega(i);
        
        % Sum all torques and calculate acceleration
        I = 1/12 * L^2; % Simple moment of inertia for a rod
        alpha(i) = (tau_magnetic + tau_spring + tau_damping) / I;
    end
    
    dY = [omega; alpha];
end

function T = compute_transform_matrices(theta, link_length)
    % Compute 2D transformation matrices for each link
    num_links = length(theta);
    T = cell(1, num_links);
    x = 0;
    y = 0;
    current_theta = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        dx = link_length * cos(current_theta);
        dy = link_length * sin(current_theta);

        % Transformation matrix (2D rotation and translation)
        T{i} = [cos(current_theta), -sin(current_theta), x; 
                sin(current_theta),  cos(current_theta), y; 
                0,                  0,                  1];
        % Update position
        x = x + dx;
        y = y + dy;
    end
end