classdef cost_function

    properties

        RS = RobotState;

    end

    methods

        function cost = moment_equilibrium(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            psi = x(1:num_links);
            theta_M = x(num_links+1:2*num_links);
            theta_link = x(2*num_links+1:end);
            r_ext = [r; -0.00165];
            M = psi * link_length * cross_section_area;

            link_center = obj.RS.set_link_center(num_links, link_length, theta_link);

            mag_to_link = link_center - repmat(r_ext, 1, num_links);


            theta = 0;
            theta_M_each = zeros(1, num_links);

            eq1 = zeros(1, num_links);
            eq2 = zeros(1, num_links);
            tm = zeros(1, num_links);
            ts = zeros(1, num_links);
            tsj = zeros(1, num_links);
            cost = 0;
            for k = 1:num_links

                theta = theta + theta_link(k);
                theta_M_each(k) = theta + theta_M(k);
                m_vec = M(k) * [cos(theta_M_each(k)); sin(theta_M_each(k)); 0];
                tau_magnet = norm(cross(m_vec, EM2.Cal_B_Field(mag_to_link(:,k))));

                tau_spring = 0;

                if k > 1
                    tau_spring = tau_spring - k_spring(k-1)*(theta_link(k) - theta_link(k-1));
                end

                if k < num_links
                    tau_spring = tau_spring + k_spring(k) * (theta_link(k+1) - theta_link(k));
                end

                tm(k) = tau_magnet;
                ts(k) = tau_spring;

                eq1(k) = tau_magnet + tau_spring;

                %                 cost = cost + eq(k)^2;

            end

            tau_magnet_init = zeros(1, num_links);

            for i = 1:num_links

                m_vec_init = M(i) * [cos(theta_M(i)); sin(theta_M(i)); 0];
                r_init = [link_length/2 + link_length * (i-1); 0] - r_ext;
                tau_magnet_init(i) = norm(cross(m_vec_init, EM2.Cal_B_Field(r_init)));

            end

            tsj(1) = k_spring(1) * theta_link(1);
            eq2(1) = tsj(1) + tm(1);

            for k = 1:num_links-1

                tsj(k+1) = k_spring(k+1) * (theta_link(k+1) - theta_link(k));
                eq2(k+1) = tsj(k+1) + tm(k+1) + tm(k);

            end

            eq3 = zeros(1, num_links);
            eq3(1) = tau_magnet_init(1);

            for k = 1:num_links-1

                eq3(k+1) = tau_magnet_init(k) + tau_magnet_init(k+1);

            end

            cost = norm((eq2 - eq3).^2);
            cost = cost *1e+2;

        end

        function [c, ceq] = nonlcon(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            c = [];
            ceq = [];

            theta_link = x(2*num_links+1:end);

            joint_pos = obj.RS.set_joint_pos(num_links, link_length, theta_link);
            x_joint_pos = joint_pos(1, :);
            %             y_joint_pos = joint_pos(2, :);

            ceq = joint_pos(2,end); % y coordinate of the end joint is fixed on y-axis

            c_rest = zeros(num_links-1, 1);

            for i = 1:num_links-1

                c_rest(i) = x_joint_pos(i) - x_joint_pos(i+1);

            end

            c2 = - joint_pos(2, end) + 0.003;

            %             c = [c_rest; c2];
            c = c2;

        end


        function tau = get_tau(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            psi = x(1:num_links);
            theta_M = x(num_links+1:2*num_links);
            theta_link = x(2*num_links+1:end);
            r_ext = [r; -0.00165];
            M = psi * link_length * cross_section_area;

            link_center = obj.RS.set_link_center(num_links, link_length, theta_link);

            mag_to_link = link_center - repmat(r_ext, 1, num_links);


            theta = 0;
            theta_M_each = zeros(1, num_links);

            eq1 = zeros(1, num_links);
            eq2 = zeros(1, num_links);
            tm = zeros(1, num_links);
            ts = zeros(1, num_links);
            tsj = zeros(1, num_links);

            cost = 0;


            for k = 1:num_links

                theta = theta + theta_link(k);
                theta_M_each(k) = theta + theta_M(k);
                m_vec = M(k) * [cos(theta_M_each(k)); sin(theta_M_each(k)); 0];
                tau_magnet = norm(cross(m_vec, EM2.Cal_B_Field(mag_to_link(:,k))));

                tau_spring = 0;

                if k > 1
                    tau_spring = tau_spring - k_spring(k-1)*(theta_link(k) - theta_link(k-1));
                end

                if k < num_links
                    tau_spring = tau_spring + k_spring(k) * (theta_link(k+1) - theta_link(k));
                end

                tm(k) = tau_magnet;
                ts(k) = tau_spring;

                eq1(k) = tau_magnet + tau_spring;

                %                 cost = cost + eq(k)^2;

            end

            tau_magnet_init = zeros(1, num_links);

            for i = 1:num_links

                m_vec_init = M(i) * [cos(theta_M(i)); sin(theta_M(i)); 0];
                r_init = [link_length/2 + link_length * (i-1); 0] - r_ext;
                tau_magnet_init(i) = norm(cross(m_vec_init, EM2.Cal_B_Field(r_init)));

            end

            tsj(1) = k_spring(1) * theta_link(1);
            eq2(1) = tsj(1) + tm(1);

            for k = 1:num_links-1

                tsj(k+1) = k_spring(k+1) * (theta_link(k+1) - theta_link(k));
                eq2(k+1) = tsj(k+1) + tm(k+1) + tm(k);

            end

            eq3 = zeros(1, num_links);
            eq3(1) = tau_magnet_init(1);

            for k = 1:num_links-1

                eq3(k+1) = tau_magnet_init(k) + tau_magnet_init(k+1);

            end

            cost = norm((eq2 - eq3).^2);

            tau = [tau_magnet_init; eq2; eq3; norm(eq2), norm(eq3), cost, zeros(1, num_links-3)];

        end

        %-----------------------------another cost function--------------------------------%

        function cost = Max_Bending(obj, x, num_links, link_length, cross_section_area, r, k_spring, EM2)

            psi_init = x(1:num_links);
            theta_M = x(num_links+1:end);

            theta = RobotState;

        end

    end

end