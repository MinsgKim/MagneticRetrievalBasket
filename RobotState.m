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

        function theta_link = Get_Link_Angle(~, num_links, link_length, M, theta_M, r, k_spring, EM2)

            % requisites for calculating torques
            theta_link = zeros(1, num_links);
            joint_pos = zeros(2, num_links+1);
            link_center = zeros(2, num_links);
            theta_accum = 0; % required for calculating a center of each link

            r_ext = [r; -0.00165];

            % data storage
            T_m = zeros(1, num_links);
            T_s = zeros(1, num_links);
            T_sum = zeros(1, num_links);

            % check if it is converged or not

            maxIter = 1e+3;
            isConverged  = false(1, num_links);
            torqueTolerance = 1e-8;

            dtheta = 0.01; % about 0.0573 degree

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
%                     disp(j)
%                     disp(T_m_eq)
%                     disp(T_s(j))
                    T_sum(j) = T_m_eq - T_s(j);
%                     disp(T_sum(j))

                end

                theta_accum = 0;

                for j = 1:num_links
                    % 이미 수렴하여 멈춘 링크라면 넘어감
                    if isConverged(j)
                        continue;
                    end

                    % netTorque가 양수 -> +dtheta
                    if T_sum(j) > 0
                        theta_link(j) = theta_link(j) + dtheta;

                        % 만약 이전 사이클에서 음수였다가 이번에 양수가 됐다면(부호 변화)
                        % 혹은 토크가 0 근처라면, 이제 멈춤으로 처리
                        % (부호 바뀌는 판정은 필요에 따라 세부적으로 구현)
                        %
                        % 예시: 절댓값이 매우 작으면 멈춤
                        if abs(T_sum(j)) < torqueTolerance
                            isConverged(j) = true;
                        end

                        % netTorque가 음수 -> -dtheta
                    elseif T_sum(j) < 0
                        theta_link(j) = theta_link(j) - dtheta;

                        % 부호 변화 혹은 토크가 0 근처라면 멈춤
                        if abs(T_sum(j)) < torqueTolerance
                            isConverged(j) = true;
                        end

                    else
                        % 토크가 정확히 0 근처라면(또는 0이 되면) 이미 평형점
                        isConverged(j) = true;
                    end
                end

                if all(isConverged)
                    fprintf("모든 링크가 수렴하여 회전을 멈췄습니다. (iter = %d)\n", iter);
                    break;
                end

                if i == maxIter
                    warning("최대 반복 횟수에 도달하였습니다. 수렴하지 않았을 수 있습니다.");
                end
            end

        end

        function tau = Get_Tau(~, num_links, link_length, M, theta_M, r, k_spring, EM2)

                        % requisites for calculating torques
            theta_link = zeros(1, num_links);
            joint_pos = zeros(2, num_links+1);
            link_center = zeros(2, num_links);
            theta_accum = 0; % required for calculating a center of each link

            r_ext = [r; -0.00165];

            % data storage
            T_m = zeros(1, num_links);
            T_s = zeros(1, num_links);
            T_sum = zeros(1, num_links);

            % check if it is converged or not

            maxIter = 1e+6;
            isConverged  = false(1, num_links);
            torqueTolerance = 1e-10;

            dtheta = 0.001; % about 0.0573 degree

            for i=1:maxIter

                for j=1:num_links

                    % setup somethings
                    theta_accum = theta_accum + theta_link(j);
                    joint_pos(:,j+1) = joint_pos(:,j) + link_length * [cos(theta_accum); sin(theta_accum)];
                    link_center(:,j) = (joint_pos(:,j+1)+joint_pos(:,j))/2;

                    % calculate magnetic torque first:
                    M_vec = M(j) * [cos(theta_accum(j)+theta_M(j)); sin(theta_accum(j)+theta_M(j)); 0];
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
                    T_sum(j) = T_m_eq + T_s(j);

                end

                theta_accum = 0;

                for j = 1:num_links
                    % 이미 수렴하여 멈춘 링크라면 넘어감
                    if isConverged(j)
                        continue;
                    end

                    % netTorque가 양수 -> +dtheta
                    if T_sum(j) > 0
                        theta_link(j) = theta_link(j) + dtheta;

                        % 만약 이전 사이클에서 음수였다가 이번에 양수가 됐다면(부호 변화)
                        % 혹은 토크가 0 근처라면, 이제 멈춤으로 처리
                        % (부호 바뀌는 판정은 필요에 따라 세부적으로 구현)
                        %
                        % 예시: 절댓값이 매우 작으면 멈춤
                        if abs(T_sum(j)) < torqueTolerance
                            isConverged(j) = true;
                        end

                        % netTorque가 음수 -> -dtheta
                    elseif T_sum(j) < 0
                        theta_link(j) = theta_link(i) - dtheta;

                        % 부호 변화 혹은 토크가 0 근처라면 멈춤
                        if abs(T_sum(j)) < torqueTolerance
                            isConverged(j) = true;
                        end

                    else
                        % 토크가 정확히 0 근처라면(또는 0이 되면) 이미 평형점
                        isConverged(j) = true;
                    end
                end

                if all(isConverged)
                    fprintf("모든 링크가 수렴하여 회전을 멈췄습니다. (iter = %d)\n", iter);
                    break;
                end

                if iter == maxIter
                    warning("최대 반복 횟수에 도달하였습니다. 수렴하지 않았을 수 있습니다.");
                end
            end

            tau = [T_m; T_s; T_sum];


        end

    end

end