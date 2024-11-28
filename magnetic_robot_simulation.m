% Main script for magnetic robot simulation
function Y = magnetic_robot_simulation(x)

    x(1) = M;
    x(2) = theta_M;
    x(3) = B;
    x(4) = link_length;
    % Parameters
%     num_links = 3;
%     link_length = 0.01; % 1cm links
%     k_spring = [0.01, 0.01]; % Spring constants between links [N*m/rad]
%     M = [1, 1, 1]; % Magnetization magnitude for each link [A*m^2]
%     theta_M = [0, pi/2, 0]; % Initial magnetization angles relative to link orientation [rad]

    % rotational spring constant of PDMS
    E_pdms = 2e06;  % Young's modulus
    v_pdms = 0.5;   % Poisson's ratio
    mu_pdms = 0.67e06; % shear elastic modulus (Neo-Hookean model)
    
    % Initial configuration
    theta_init = zeros(1, num_links); % Initial link angles [rad]
    
    % Time parameters
    t_span = [0 10];
    dt = 0.01;
    t = t_span(1):dt:t_span(2);
    
    % External magnetic field parameters
%     B_magnitude = 0.01; % Tesla
%     B_angle = pi/4; % rad
%     B = B_magnitude * [cos(B_angle); sin(B_angle)];
    
    % Solve equations of motion using ode45
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, Y] = ode45(@(t,y) robot_dynamics(t, y, k_spring, x(1), x(2), x(3), x(4)), ...
                    t_span, [theta_init, zeros(1,num_links)], options);
    

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