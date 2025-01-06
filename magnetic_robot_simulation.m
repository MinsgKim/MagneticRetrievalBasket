classdef magnetic_robot_simulation

    properties
        RK = Robot_Kinematics;
        EM = External_Magnet;
    end

    methods

        function cost = objective_function(obj, x, num_links)
            % cost function: maximize summation of middle links -> bending a lot
            % optim variables
            M = x(1:num_links);
            theta_M = x(num_links+1:2*num_links);
            r = x(2*num_links+1);
            link_length = x(end);
            cross_section_area = 0.0033 * 0.0005; % cross sectional area
            M = M * link_length * cross_section_area;

            gamma = 2.0;

            % magnetic robot simulation
            [T_actual, ~] = obj.RK.simulate_robot_transform(num_links, M, theta_M, r, link_length, obj.EM);

            % calculate the sum of middle links
            x_positions = zeros(num_links-2, 1);
            y_positions = zeros(num_links-2, 1);
            for i = 2:num_links-1
                x_positions(i-1) = T_actual{i}(1, 3);
                y_positions(i-1) = T_actual{i}(2, 3);
            end
            %             cost = sum(x_positions);
            x_mid = x_positions(3:4);
            y_mid = y_positions(3:4);
            x_sum = sum(x_mid);
            y_sum = sum(y_mid);

            cost = x_sum + gamma * y_sum;

            % To maximize, mimimize the cost
            cost = -cost;
        end

        function [c, ceq] = nonlcon_position_constraints(obj, x, num_links)
            c = [];
            ceq = [];
            % nonlinear constraints setup
            % optim variables
            M = x(1:num_links);
            theta_M = x(num_links+1:2*num_links);
            r = x(2*num_links+1);
            link_length = x(end);
            cross_section_area = 0.0033 * 0.0005; % cross sectional area
            M = M * link_length * cross_section_area;

            % magnetic robot simulation
            [T_actual, ~] = obj.RK.simulate_robot_transform(num_links, M, theta_M, r, link_length, obj.EM);

            % extract the position of each link
            positions = zeros(2, num_links);
            for i = 1:num_links
                positions(:, i) = T_actual{i}(1:2, 3);
            end

            % equation constraints: x position of 1st and last links are same
            %             ceq_1 = positions(1, 1);
            %             ceq_2 = positions(2, 1);
            %
            %             ceq = [ceq_1; ceq_2];
            % nonequation constraints:
            %             c = positions(1, 1) - positions(1, 2:end-1) + 1e-6; % minute tolerance 1e-6
            %             link_distances = sqrt(sum(diff(positions, 1, 2).^2));
            %             c = link_distances - link_length; % 각 링크 간 거리가 link_length와 일치하도록 강제
            %             y_con = 0.002 - positions(2, end); % y coord. of last link > 2 mm (not negative)
            %             c = [c(:); y_con];
            x_con1 = positions(1, end) - 0.003;
            %             x_con2 = positions(1, end) + 0.002;
            %             c = [x_con1, x_con2];
            c = x_con1;
        end

%--------------------quasi-static equilibrium equation----------------------

        function cost = objective_static(obj, x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)
            % 전역변수로 best_fval, best_theta 선언
            global best_fval best_theta theta_test2

            % 1) x = [psi_1..psi_n, thetaM_1..thetaM_n]
            gamma = 2.0; % y 항에 대한 가중치
            psi = x(1:num_links);
            thetaM = x(num_links+1 : 2*num_links);

            % 2) A/m -> A·m^2
            M = psi .* (link_length * cross_section_area);

            % 3) 정적 해석 -> theta_eq
            theta_eq = obj.RK.solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring);

            % 4) 중간 링크 x좌표 합
            pos = obj.RK.compute_link_positions2(theta_eq, link_length);
            %     if num_links>2
            %         x_mid = pos(1, 2:end-1);
            %         y_mid = pos(2, 2:end-1);
            %     else
            %         x_mid = 0;
            %         y_mid = 0;
            %     end
            %     x_sum = sum(x_mid);
            %     y_sum = sum(y_mid);
            %
            %     y_end = pos(2,end);

            x_mid = pos(1, 3:4);
            y_mid = pos(2, 3:4);
            x_sum = sum(x_mid);
            y_sum = sum(y_mid);

            cost = -(x_sum + gamma * y_sum);
            
            theta_test2 = theta_eq;

            % => cost를 최소화 => x_sum + gamma*y_end 최대화
            % 5) "현재 cost가 더 좋으면" -> best_fval 갱신, best_theta 갱신
            if cost < best_fval
                best_fval = cost;
                best_theta = theta_eq;
            end

        end

        function [c, ceq] = constraint_static(obj, x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)
            c=[];
            ceq=[];

            psi = x(1:num_links);
            thetaM = x(num_links+1 : 2*num_links);

            % A/m -> A·m^2
            M = psi .* (link_length * cross_section_area);

            % 정적 해석
            theta_eq = obj.RK.solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring);

            % 링크 끝단 위치
            pos = obj.RK.compute_link_positions2(theta_eq, link_length);

            % --- (1) 첫 링크를 (0,0)에 고정 ---
            ceq_first_x = pos(1,1);  % = 0
            ceq_first_y = pos(2,1);  % = 0

            % --- (2) 링크 간 거리(부등호) ---
            link_distances = sqrt(sum(diff(pos,1,2).^2,1));
            c_link_dist = link_distances - link_length;  % <= 0

            % --- (3) 마지막 링크 x in [-1 mm, +1 mm], y >= 3 mm (예시) ---
            x_end = pos(1,end);
            y_end = pos(2,end);

            x_con1 = x_end - 0.003;   % <=0  -> x_end <= +1 mm
            x_con2 = -x_end + 0.01;  % <=0  -> x_end >= -1 mm
            %     y_con  = 0.003 - y_end;   % <=0  -> y_end >= 3 mm

            % 등가 제약(모두 0이어야 함)
            %     ceq = [ceq_first_x; ceq_first_y];


            % 부등호 제약(모두 <= 0)
            c = [x_con1];
        end



    end

end