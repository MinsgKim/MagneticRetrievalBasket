classdef cost_function

    properties

        RS = RobotState;

    end

    methods

        function cost = Max_Bending(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            psi_init = x(1:num_links);
            theta_M = x(num_links+1:end);
            M = psi_init * link_length * cross_section_area;

            theta = obj.RS.Get_Link_Angle(num_links, link_length, M, theta_M, r, k_spring, EM2);

            y_joint_pos = zeros(1, num_links + 1);
            theta_sum = 0;

            for i=1:num_links
                
                theta_sum = theta_sum + theta(i);
                y_joint_pos(i+1) = y_joint_pos(i) + link_length * sin(theta_sum);

            end
            
            if num_links/2 == 0
                
                cost = y_joint_pos(num_links/2+1);

            else
                
                cost = (y_joint_pos((num_links+1)/2) + y_joint_pos((num_links+1)/2+1))/2;

            end

            cost = -cost;

        end

        function [c, ceq] = nonlcon_ver2(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            c = [];

            psi_init = x(1:num_links);
            theta_M = x(num_links+1:end);
            M = psi_init * link_length * cross_section_area;

            theta = obj.RS.Get_Link_Angle(num_links, link_length, M, theta_M, r, k_spring, EM2);

            y_joint_pos = zeros(1, num_links + 1);
            theta_sum = 0;

            for i=1:num_links
                
                theta_sum = theta_sum + theta(i);
                y_joint_pos(i+1) = y_joint_pos(i) + link_length * sin(theta_sum);

            end

            ceq = y_joint_pos(end);
%             c = y_joint_pos(5) - 0.004;

        end

    end

end