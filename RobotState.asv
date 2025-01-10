classdef RobotState

    properties

    end

    methods

        function joint_pos = set_joint_pos(~, num_links, link_length, theta_link)

            theta = 0;
            x_joint_pos = zeros(1, num_links+1);
            y_joint_pos = zeros(1, num_links+1);

            for i = 1:num_links

                theta = theta + theta_link(i);
                x_joint_pos(i+1) = x_joint_pos(i) + link_length * cos(theta);
                y_joint_pos(i+1) = y_joint_pos(i) + link_length * sin(theta);

            end

            joint_pos = [x_joint_pos; y_joint_pos];

        end

        function link_center = set_link_center(~, num_links, link_length, theta_link)

            x_joint_pos = zeros(1, num_links+1);
            y_joint_pos = zeros(1, num_links+1);
            x_link_center = zeros(1, num_links);
            y_link_center = zeros(1, num_links);

            theta = 0;

            for i = 1:num_links

                theta = theta + theta_link(i);
                x_joint_pos(i+1) = x_joint_pos(i) + link_length * cos(theta);
                y_joint_pos(i+1) = y_joint_pos(i) + link_length * sin(theta);
                x_link_center(i) = (x_joint_pos(i+1) + x_joint_pos(i))/2;
                y_link_center(i) = (y_joint_pos(i+1) + y_joint_pos(i))/2;

            end

            link_center = [x_link_center; y_link_center];

        end

        function draw_plot(obj, num_links, link_length, theta_link, theta_M)

            figure;
            hold on

            joint_pos = obj.set_joint_pos(num_links, link_length, theta_link);
            link_center = obj.set_link_center(num_links, link_length, theta_link);

            plot(joint_pos(1,:), joint_pos(2,:), '-o', 'LineWidth', 2);
            xlabel('X position')
            xlabel('X position')
            grid on
            axis equal

            th = 0;
            for i = 1:num_links

                th = th + theta_link(i);
                theta2 = th + theta_M(i);
                x_magntz = cos(theta2);
                y_magntz = sin(theta2);

                quiver(link_center(1,i), link_center(2,i), x_magntz, y_magntz, link_length, 'Color', 'r', 'LineWidth', 1.25)

            end

        end

        %------------------------new function----------------------------%

        function theta_link = Get_Link_Angle(~, num_links, link_length, M, r, k_spring, EM2)

            

        end

        function tau = Get_Tau(~, num_links, link_length, M, r, k_spring, EM2)

            

        end

    end

end