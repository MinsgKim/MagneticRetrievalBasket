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
        %
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

        function theta_link = Get_Link_Angle(~, num_links, link_length, M, theta_M, r, k_spring, EM2)

            % requisites for calculating torques
            theta_link = zeros(1, num_links);
            joint_pos = zeros(2, num_links+1);
            link_center = zeros(2, num_links);
            theta_accum = 0; % required for calculating a center of each link

            r_ext = [r; -0.0015];

            % data storage
            T_m = zeros(1, num_links);
            T_s = zeros(1, num_links);
            T_sum = zeros(1, num_links);

            % check if it is converged or not

            maxIter = 2e+3;
            isConverged  = false(maxIter, num_links);
            %             torqueTolerance = 1e-10;

            T_m_add_storage = zeros(maxIter, num_links);
            T_s_storage = zeros(maxIter, num_links);
            T_sum_storage = zeros(maxIter, num_links);

            dtheta = 0.002; % about 0.0573 degree

            for i=1:maxIter

                for j=1:num_links

                    % setup somethings
                    theta_accum = theta_accum + theta_link(j);
                    joint_pos(:,j+1) = joint_pos(:,j) + link_length * [cos(theta_accum); sin(theta_accum)];
                    link_center(:,j) = (joint_pos(:,j+1)+joint_pos(:,j))/2;

                    % calculate magnetic torque first:
                    M_vec = M(j) * [cos(theta_accum+theta_M(j)); sin(theta_accum+theta_M(j)); 0];
                    r_vec = link_center(:,j) - r_ext;
                    B_Field = EM2.Cal_B_Field(r_vec);
                    tau_mag = cross(M_vec, B_Field);
                    T_m(j) = tau_mag(3);

                    if j > 1
                        T_m_eq = T_m(j) - T_m(j-1);
                        T_s(j) = k_spring(j)*(theta_link(j)-theta_link(j-1));

                    else
                        T_m_eq = T_m(j);
                        T_s(j) = k_spring(j)*theta_link(j);

                    end
                    % calculate net Torque
                    T_sum(j) = T_m_eq - T_s(j);

                    T_m_add_storage(i,j) = T_m_eq;
                    T_s_storage(i,j) = T_s(j);
                    T_sum_storage(i,j) = T_sum(j);

                end

                theta_accum = 0;

                if i == 1

                    if isConverged(i,1) == false
                        if T_sum_storage(i,1) > 0
                            theta_link(1) = theta_link(1) + dtheta;
                        elseif T_sum_storage(i,1) < 0
                            theta_link(1) = theta_link(1) - dtheta;
                        else
                            isConverged(i,1) = true;
                        end
                    end

                else
                    if isConverged(i-1,1) == false
                        if T_sum_storage(i,1) > 0
                            if T_sum_storage(i-1, 1) <= 0
                                isConverged(i,1) = true;
                            else
                                theta_link(1) = theta_link(1) + dtheta;
                            end

                        elseif T_sum_storage(i,1) < 0
                            if T_sum_storage(i-1,1) >= 0
                                isConverged(i,1) = true;
                            else
                                theta_link(1) = theta_link(1) - dtheta;
                            end

                        else
                            isConverged(i,1) = true;

                        end

                    elseif isConverged(i-1,1) == true
                        theta_link(1) = theta_link(1);
                        isConverged(i,1) = true;
                    end

                end

                for j=2:num_links

                    if i == 1
                        if isConverged(i,j) == false
                            if T_sum_storage(i,j) > 0
                                theta_link(j) = theta_link(j) + dtheta;

                            elseif T_sum_storage(i,j) < 0
                                theta_link(j) = theta_link(j) - dtheta;

                            else
                                isConverged(i,j) = true;
                            end

                        else
                            theta_link(j) = theta_link(j);
                            isConverged(i,j) = true;
                        end

                    else
                        if isConverged(i-1,j) == false
                            if T_sum_storage(i,j) > 0
                                if T_sum_storage(i-1, j) <= 0
                                    if isConverged(i,j-1) == true
                                        isConverged(i,j) = true;
                                    else
                                        theta_link(j) = theta_link(j) - dtheta;
                                    end
                                else
                                    theta_link(j) = theta_link(j) + dtheta;

                                end

                            elseif T_sum_storage(i,j) < 0
                                if T_sum_storage(i-1,j) >= 0
                                    if isConverged(i,j-1) == true
                                        isConverged(i,j) = true;
                                    else
                                        theta_link(j) = theta_link(j) + dtheta;
                                    end
                                else
                                    theta_link(j) = theta_link(j) - dtheta;
                                end

                            else
                                isConverged(i,j) = true;

                            end

                        elseif isConverged(i-1,j) == true
                            theta_link(j) = theta_link(j);
                            isConverged(i,j) = true;
                        end

                    end

                end

                if all(isConverged(i,:))
                    %                     fprintf("All links are converged. Stop rotation. (iter = %d)\n", i);
                    break;
                end

                if i == maxIter
                    warning("Reached maximum iteration. It might not be converged.");
                end
            end

        end

        function [T_m, T_s, T_sum, isConverged] = Get_Tau(~, num_links, link_length, M, theta_M, r, k_spring, EM2)

            % requisites for calculating torques
            theta_link = zeros(1, num_links);
            joint_pos = zeros(2, num_links+1);
            link_center = zeros(2, num_links);
            theta_accum = 0; % required for calculating a center of each link

            r_ext = [r; -0.0015];

            % data storage
            T_m = zeros(1, num_links);
            T_s = zeros(1, num_links);
            T_sum = zeros(1, num_links);

            % check if it is converged or not

            maxIter = 2e+3;
            isConverged  = false(maxIter, num_links);
            %             torqueTolerance = 1e-10;

            T_m_add_storage = zeros(maxIter, num_links);
            T_s_storage = zeros(maxIter, num_links);
            T_sum_storage = zeros(maxIter, num_links);

            dtheta = 0.002; % about 0.0573 degree

            for i=1:maxIter

                for j=1:num_links

                    % setup somethings
                    theta_accum = theta_accum + theta_link(j);
                    joint_pos(:,j+1) = joint_pos(:,j) + link_length * [cos(theta_accum); sin(theta_accum)];
                    link_center(:,j) = (joint_pos(:,j+1)+joint_pos(:,j))/2;

                    % calculate magnetic torque first:
                    M_vec = M(j) * [cos(theta_accum+theta_M(j)); sin(theta_accum+theta_M(j)); 0];
                    r_vec = link_center(:,j) - r_ext;
                    B_Field = EM2.Cal_B_Field(r_vec);
                    tau_mag = cross(M_vec, B_Field);
                    T_m(j) = tau_mag(3);

                    if j > 1
                        T_m_eq = T_m(j) - T_m(j-1);
                        T_s(j) = k_spring(j)*(theta_link(j)-theta_link(j-1));

                    else
                        T_m_eq = T_m(j);
                        T_s(j) = k_spring(j)*theta_link(j);

                    end
                    % calculate net Torque
                    T_sum(j) = T_m_eq - T_s(j);

                    T_m_add_storage(i,j) = T_m_eq;
                    T_s_storage(i,j) = T_s(j);
                    T_sum_storage(i,j) = T_sum(j);

                end

                theta_accum = 0;

                if i == 1

                    if isConverged(i,1) == false
                        if T_sum_storage(i,1) > 0
                            theta_link(1) = theta_link(1) + dtheta;
                        elseif T_sum_storage(i,1) < 0
                            theta_link(1) = theta_link(1) - dtheta;
                        else
                            isConverged(i,1) = true;
                        end
                    end

                else
                    if isConverged(i-1,1) == false
                        if T_sum_storage(i,1) > 0
                            if T_sum_storage(i-1, 1) <= 0
                                isConverged(i,1) = true;
                            else
                                theta_link(1) = theta_link(1) + dtheta;
                            end

                        elseif T_sum_storage(i,1) < 0
                            if T_sum_storage(i-1,1) >= 0
                                isConverged(i,1) = true;
                            else
                                theta_link(1) = theta_link(1) - dtheta;
                            end

                        else
                            isConverged(i,1) = true;

                        end

                    elseif isConverged(i-1,1) == true
                        theta_link(1) = theta_link(1);
                        isConverged(i,1) = true;
                    end

                end

                for j=2:num_links

                    if i == 1
                        if isConverged(i,j) == false
                            if T_sum_storage(i,j) > 0
                                theta_link(j) = theta_link(j) + dtheta;

                            elseif T_sum_storage(i,j) < 0
                                theta_link(j) = theta_link(j) - dtheta;

                            else
                                isConverged(i,j) = true;
                            end

                        else
                            theta_link(j) = theta_link(j);
                            isConverged(i,j) = true;
                        end

                    else
                        if isConverged(i-1,j) == false
                            if T_sum_storage(i,j) > 0
                                if T_sum_storage(i-1, j) <= 0
                                    if isConverged(i,j-1) == true
                                        isConverged(i,j) = true;
                                    else
                                        theta_link(j) = theta_link(j) - dtheta;
                                    end
                                else
                                    theta_link(j) = theta_link(j) + dtheta;

                                end

                            elseif T_sum_storage(i,j) < 0
                                if T_sum_storage(i-1,j) >= 0
                                    if isConverged(i,j-1) == true
                                        isConverged(i,j) = true;
                                    else
                                        theta_link(j) = theta_link(j) + dtheta;
                                    end
                                else
                                    theta_link(j) = theta_link(j) - dtheta;
                                end

                            else
                                isConverged(i,j) = true;

                            end

                        elseif isConverged(i-1,j) == true
                            theta_link(j) = theta_link(j);
                            isConverged(i,j) = true;
                        end

                    end

                end

                if all(isConverged(i,:))
                    %                     fprintf("All links are converged. Stop rotation. (iter = %d)\n", i);
                    break;
                end

                if i == maxIter
                    warning("Reached maximum iteration. It might not be converged.");
                end
            end

            T_m = T_m_add_storage;
            T_s = T_s_storage;
            T_sum = T_sum_storage;

        end

    end

end