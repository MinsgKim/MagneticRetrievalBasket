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
            cost = abs(x_positions(3));

            % To maximize, mimimize the cost
            cost = -cost;
        end

        function [c, ceq] = nonlcon_position_constraints(obj, x, num_links)
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
            ceq_1 = positions(1, 1);
            ceq_2 = positions(2, 1);
            
            ceq = [ceq_1; ceq_2];
            % nonequation constraints:
            %             c = positions(1, 1) - positions(1, 2:end-1) + 1e-6; % minute tolerance 1e-6
%             link_distances = sqrt(sum(diff(positions, 1, 2).^2));
%             c = link_distances - link_length; % 각 링크 간 거리가 link_length와 일치하도록 강제
%             y_con = 0.002 - positions(2, end); % y coord. of last link > 2 mm (not negative)
%             c = [c(:); y_con];
            x_con1 = positions(1, end) - 0.001;
            x_con2 = positions(1, end) + 0.001;
            c = [x_con1, x_con2];
        end


        function cost = objective_static(obj, x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)

        end

        function 


    end

end