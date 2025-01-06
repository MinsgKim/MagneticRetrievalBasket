classdef Robot_Kinematics

    properties

        k_spring = [];         % spring coefficient of joint (PDMS or Ecoflex)
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
            t_span = [0, 10]; % simulation time
            options_ode = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', 0.1);

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
            Kspring = 3e-5 .* ones(1, num_links_ - 1);
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

        function plot_robot(~, T_actual, theta, theta_M, link_length, r_init)
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

            th = 0;
            for i = 1:num_links__

                th = th + theta(i);
                theta2 = th + theta_M(i);
                x_center = (x_positions(i+1)+x_positions(i))/2;
                y_center = (y_positions(i+1)+y_positions(i))/2;
                x_magntz = sin(theta2);
                y_magntz = cos(theta2);

                quiver(x_center, y_center, x_magntz, y_magntz, link_length,'Color', 'r')

            end

            hold off;
        end



        function theta_eq = solve_static_equilibrium(obj, num_links, M, theta_M, r, link_length, EM, k_spring)
            %----------------------------------------------------------------------
            % 정적(Quasi-static) 균형 해(theta_eq)를 찾는 함수
            % 각 링크에 걸리는 순토크를 0으로 만드는 각도 벡터 theta를 fsolve로 구한다.
            %
            % INPUT
            %   num_links   : 링크 개수
            %   M           : 각 링크의 자화 세기(스칼라)
            %   theta_M     : 각 링크의 자화 방향(라디안)
            %   r           : 외부 자석과 로봇 끝단 사이 거리
            %   link_length : 링크 길이
            %   EM          : External_Magnet 객체
            %
            % OUTPUT
            %   theta_eq    : (1 x num_links) 정적 균형에서의 링크 각도
            %----------------------------------------------------------------------

            % 1) obj 프로퍼티 세팅
            obj.num_links = num_links;
            obj.M = M;
            obj.theta_M = theta_M;
            obj.r = r;
            obj.link_length = link_length;
            obj.EM = EM;
            obj.k_spring = k_spring;

            % 초기 추정값 (적당히 0 혹은 작은 난수 등)
            rng(0)
            theta_init = (pi/4) * randn(1, num_links);

            % fsolve 옵션 설정
            options = optimoptions('fsolve',...
                'Display','None',...       % 중간 과정 표시
                'MaxIterations',1000,...
                'MaxFunctionEvaluations',1e5,...
                'FunctionTolerance',1e-12);

            % fsolve 실행
            theta_eq = fsolve(@(theta) obj.equilibrium_equations(theta), ...
                theta_init, options);

        end

        function eq = equilibrium_equations(obj, theta)
            num_links_ = obj.num_links;
            M_        = obj.M;
            theta_M_  = obj.theta_M;
            r_        = obj.r;
            link_length_ = obj.link_length;
            EM_       = obj.EM;

            %-----------------------------------------------------------
            % 입력 theta에 대해, 각 링크에 작용하는
            %   (자화 토크 + 스프링 토크)의 합이 0이어야 함
            %   -> eq(i) = 0
            %-----------------------------------------------------------
            eq = zeros(1, num_links_);


            % 스프링 계수(예: 링크 개수-1개만큼)
            Kspring = obj.k_spring;

            % 링크들의 2D 위치 계산
            positions = obj.compute_link_positions2(theta, link_length_);

            % 각 링크에 대한 순토크
            for i = 1:num_links_

                %----------------------------
                % 1) 자기장 계산
                %----------------------------
                % 외부 자석 기준으로 r_vec 설정
                r_vec = positions(:, i) - [-0.00165; num_links_*link_length_ + r_];
                % B_ext: External_Magnet 클래스의 Dipole 식으로 계산
                B_ext = EM_.Cal_B(r_vec);

                % 링크 자화 모멘트 m_i
                theta_m = theta(i) + theta_M_(i);  % 자화 방향 = (링크 각도 + 자화 편차)
                m_i = M_(i) * [sin(theta_m); cos(theta_m); 0];

                % 2D 상에서 tau_magnetic = (m x B)의 z성분
                tau_magnetic = m_i(2)*B_ext(1) - m_i(1)*B_ext(2);

                %----------------------------
                % 2) 스프링 토크
                %   - (i-1) 링크와 (i) 링크간 상대 회전각
                %----------------------------
                tau_spring = 0;
                if i > 1
                    tau_spring = tau_spring - Kspring(i - 1)*(theta(i) - theta(i - 1));
                end
                if i < num_links_
                    tau_spring = tau_spring + Kspring(i) * (theta(i+1) - theta(i));
                end

                %----------------------------
                % (damping, inertia는 정적에서 0)
                %----------------------------
                % 최종 합 (정적이므로 = 0)
                eq(i) = tau_magnetic + tau_spring;
            end
        end

        function positions = compute_link_positions2(~, theta, link_length)
            n = length(theta);
            positions = zeros(2, n);

            x = 0; y = 0;
            current_angle = 0;
            for i = 1:n
                current_angle = current_angle + theta(i);
                x = x + link_length * sin(current_angle);
                y = y + link_length * cos(current_angle);
                positions(:, i) = [x; y];
            end
        end

    end


end
