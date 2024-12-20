classdef Robot_Kinematics

    properties

        k_spring = 9.21e-4;         % spring coefficient of joint (PDMS or Ecoflex)
        damping = 1e-6; % damping coefficient of joint
        num_links;
        M;
        theta_M;
        r;
        link_length;
        EM;


    end

    methods

        function [T_actual, theta_final] = simulate_robot_transform(obj, num_links, M, theta_M, r, link_length, EM)

            % parameter setup
            theta_init = 1e-3 * randn(1, num_links); % initial random angle of each link
            t_span = [0, 50]; % simulation time
            options_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.1);

            % Store data in object properties
            obj.num_links = num_links;
            obj.M = M;
            obj.theta_M = theta_M;
            obj.r = r;
            obj.link_length = link_length;
            obj.EM = EM;

            % Initial conditions
            Y0 = [theta_init, zeros(1, num_links)];

            % Simulation with ode15s
            [~, Y] = ode15s(@(t, y) obj.robot_dynamics(t, y), t_span, Y0, options_ode);

            % Extract the final link angle
            theta_final = Y(end, 1:num_links);

            % Calculate transformation matrix
            T_actual = obj.compute_transform_matrices(theta_final, link_length);
        end

        function dY = robot_dynamics(obj, ~, Y)
            % Access data from object properties
            num_links_ = obj.num_links;
            M_ = obj.M;
            theta_M_ = obj.theta_M;
            r_ = obj.r;
            L = obj.link_length;
            EM_ = obj.EM;

            % Initialize spring coefficient and damping
            Kspring = 9.21e-4 .* ones(1, num_links_ - 1);
            damping_coefficient = 1e-6;

            % Extract state vector
            theta = Y(1:num_links_);
            omega = Y(num_links_+1:end);

            % Initiate acceleration vector
            alpha = zeros(num_links_, 1);

            % Inertia moment (unit mass, length L)
            I = (1/12) * L^2;

            % Calculate the position of each link
            positions = obj.compute_link_positions(theta, L);

            % Calculate torque on each link
            for i = 1:num_links_
                % Vector from an external magnet to a link
                r_vec = positions(:, i) - [-0.00165; num_links_ * L + r_];

                % Calculate magnetic field on a link
                B_ext = EM_.Cal_B(r_vec);

                % Magnetic moment vector of a link
                theta_moment = theta(i) + theta_M_(i);
                m_i = M_(i) * [sin(theta_moment); cos(theta_moment); 0];

                % Calculate magnetic torque (only consider z-axis)
                tau_magnetic = m_i(2) * B_ext(1) - m_i(1) * B_ext(2);

                % Spring torque
                tau_spring = 0;
                if i > 1
                    tau_spring = tau_spring - Kspring(i - 1) * (theta(i) - theta(i - 1));
                end
                if i < num_links_
                    tau_spring = tau_spring + Kspring(i) * (theta(i + 1) - theta(i));
                end

                % Damping torque
                tau_damping = -damping_coefficient * omega(i);

                % Find angular acceleration
                alpha(i) = (tau_magnetic + tau_spring + tau_damping) / I;
            end

            % Concatenate omega and alpha for ode15s output
            dY = [omega; alpha];
        end


        function positions = compute_link_positions(obj, theta, link_length)
            % Calculate the end position of each link
%             num_links = length(theta);
            x = zeros(1, obj.num_links);
            y = zeros(1, obj.num_links);
            current_theta = 0;
            x_current = 0;
            y_current = 0;

            for i = 1:obj.num_links
                current_theta = current_theta + theta(i);
                x_current = x_current + link_length * sin(current_theta);
                y_current = y_current + link_length * cos(current_theta);
                x(i) = x_current;
                y(i) = y_current;
            end

            positions = [x; y]; % End point of each link
        end

        function T = compute_transform_matrices(~, theta, link_length)
            % calculate 2d transformation matrix
            obj.num_links = length(theta);
            T = cell(1, obj.num_links);
            x = 0;
            y = 0;
            current_theta = 0;

            for i = 1:obj.num_links
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

        function plot_robot(~, T_actual)
            figure;
            hold on;
            num_links__ = length(T_actual);
            x_positions = zeros(1, num_links__ + 1);
            y_positions = zeros(1, num_links__ + 1);
            x_positions(1) = 0;
            y_positions(1) = 0;
            for i = 1:num_links__
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

    end

end
