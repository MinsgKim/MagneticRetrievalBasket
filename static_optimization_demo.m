function static_optimization_demo()
    clc; clear; close all;

    %%% 전역 변수 선언 (최적화 중 기록용) %%%
    global best_fval best_theta
    best_fval = inf;
    best_theta = [];

    % 로봇 링크 개수
    num_links = 7;
    link_length = 2e-3; 
    cross_section_area = 0.0033 * 0.0005;
    k_spring = 3e-5 * ones(1, num_links-1);
    r_ext = 0.04;
    EM = External_Magnet();

    % ---- 최적화 변수 정의 (psi, theta_M) ----
    lb_psi = 5e3;  ub_psi = 6e4;
    lb_thetaM = -pi;  ub_thetaM = pi;
    lb = [lb_psi*ones(1,num_links), lb_thetaM*ones(1,num_links)];
    ub = [ub_psi*ones(1,num_links), ub_thetaM*ones(1,num_links)];

    % 초기값
    psi_init    = 5e3 * (7 * rand(1,num_links) + 1);
    thetaM_init = (rand(1,num_links)*2*pi - pi);
    x0 = [psi_init, thetaM_init];

    % optimized parameters storage
    num_iterations = 10;
    x_results = zeros(num_iterations, length(x0));
    cost_values = zeros(num_iterations, 1);

    % fmincon 옵션
    opts = optimoptions('fmincon',...
        'Display','iter',...
        'MaxFunctionEvaluations',1e5,...
        'MaxIterations',1000,...
        'Algorithm','sqp');

    % createOptimProblem + GlobalSearch
    problem = createOptimProblem('fmincon',...
        'x0', x0, ...
        'lb', lb, 'ub', ub, ...
        'objective', @(x) objective_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
        'nonlcon', @(x) constraint_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
        'options', opts);

    gs = GlobalSearch;

    for i=1:num_iterations
        [x_opt, fval_opt] = run(gs, problem);
        x_results(i,:) = x_opt;
        cost_values(i) = fval_opt;
        fprintf("Iteration %d completed.\n", i)
    end

    [~, best_cost] = min(cost_values);
    x_opt = x_results(best_cost, :);

    fprintf('==== Optimization Finished ====\n');
    fprintf('Best Cost (fval_opt): %.4f\n', min(cost_values));
    disp('Optimized x_opt = [psi_1..psi_n, thetaM_1..thetaM_n]:');
    disp(x_opt);

    %%% 여기서 더 이상 solve_static_equilibrium을 새로 부르지 않고,
    %%% "best_theta"를 최종 로봇 형상으로 활용
    %%% --------------------------------------------
    global best_theta
    if isempty(best_theta)
        warning('best_theta is empty! Possibly no update occurred?');
        return;
    end

    theta_final = best_theta;
    positions_opt = compute_link_positions(theta_final, link_length);

    % 검산: 중간 링크 x좌표 합
    if num_links>2
        x_mid = positions_opt(1,2:end-1);
    else
        x_mid = 0;
    end
    cost_examination = -sum(x_mid);
    fprintf('Re-check cost with best_theta = %.4f\n', cost_examination);

    %% 최종 플롯
    figure; hold on; axis equal;

    for i = 1:num_links
        r_vec = positions_opt(:,i) - [0; r_ext];
        B_local = EM.Cal_B(r_vec);
        fprintf('Link %d:  r=%.4e, |B|=%.4e\n', i, norm(r_vec), norm(B_local));
    end

    % 로봇 형상
    for i = 1:num_links-1
        plot(positions_opt(1,i:i+1), positions_opt(2,i:i+1), 'bo-','LineWidth',2 );
    end
    plot([0 positions_opt(1,1)], [0 positions_opt(2,1)], 'bo-','LineWidth',2)
    plot([0, positions_opt(1,:)], [0, positions_opt(2,:)], 'ro','MarkerSize',8,'LineWidth',2 );

    % 자화 방향 화살표
    psi_opt    = x_opt(1:num_links);
    thetaM_opt = x_opt(num_links+1 : 2*num_links);
    current_th = 0;
    for i = 1:num_links
        current_th = current_th + theta_final(i);
        th_m = theta_final(i) + thetaM_opt(i);
        if i == 1
            x_center = positions_opt(1,i)/2;
            y_center = positions_opt(2,i)/2;
        else
            x_center = (positions_opt(1,i)+positions_opt(1,i-1))/2;
            y_center = (positions_opt(2,i)+positions_opt(2,i-1))/2;
        end
        quiver( x_center, y_center, 0.5*link_length*sin(th_m), 0.5*link_length*cos(th_m), ...
            'Color',[1,0,0], 'LineWidth',1.5, 'MaxHeadSize',2 );
    end

    title('Final Robot Configuration (Best Cost)');
    xlabel('X (m)'); ylabel('Y (m)');
    grid on;
    hold off;

end


%% =========================================================================
%% (A) 목적함수: 중간 링크들의 x좌표 합을 "최대한" 크게 -> 음수 부호로 최소화
%% =========================================================================
function cost = objective_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)
    % 전역변수로 best_fval, best_theta 선언
    global best_fval best_theta

    % 1) x = [psi_1..psi_n, thetaM_1..thetaM_n]
    gamma = 2.0; % y 항에 대한 가중치
    psi = x(1:num_links);
    thetaM = x(num_links+1 : 2*num_links);

    % 2) A/m -> A·m^2
    M = psi .* (link_length * cross_section_area);

    % 3) 정적 해석 -> theta_eq
    theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring);

    % 4) 중간 링크 x좌표 합
    pos = compute_link_positions(theta_eq, link_length);
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
    % => cost를 최소화 => x_sum + gamma*y_end 최대화
    % 5) "현재 cost가 더 좋으면" -> best_fval 갱신, best_theta 갱신
    if cost < best_fval
        best_fval = cost;
        best_theta = theta_eq;
    end

end


%% =========================================================================
%% (B) 비선형 제약조건
%% =========================================================================
function [c, ceq] = constraint_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)
c=[];
ceq=[];

    psi = x(1:num_links);
    thetaM = x(num_links+1 : 2*num_links);

    % A/m -> A·m^2
    M = psi .* (link_length * cross_section_area);

    % 정적 해석
    theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring);

    % 링크 끝단 위치
    pos = compute_link_positions(theta_eq, link_length);

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

%% =========================================================================
%% (C) 정적 해석(각 링크에 토크합=0) -> fsolve
%% =========================================================================
function theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring)
    % 초기값 (약간 랜덤 or 0)
%     rng(0)
    theta_init = (pi/4) * randn(1, num_links);

    opts = optimoptions('fsolve',...
        'Display','none',...
        'MaxIterations',1000,...
        'MaxFunctionEvaluations',1e5,...
        'FunctionTolerance',1e-12);

    theta_eq = fsolve(@(theta) equilibrium_equations(theta, num_links, M, thetaM, r_ext, link_length, EM, k_spring), ...
                      theta_init, opts);

    
end

%% (C-1) 각 링크별 토크합(=0) 방정식
function eq = equilibrium_equations(theta, num_links, M, thetaM, r_ext, link_length, EM, k_spring)
    eq = zeros(1, num_links);

    % 링크 끝단 위치
    positions = compute_link_positions(theta, link_length);

    for i = 1:num_links
        % ---- 1) 외부 자기장 계산 ----
        %   pos(:,i): i번 링크 끝단 위치
        %   (예) 외부 자석이 (0, r_ext)에 있다고 단순 가정
        r_vec = positions(:, i) - [0; r_ext];
        B_ext = EM.Cal_B(r_vec);

        % ---- 2) 자기 토크 (2D에서 z성분) ----
        % m_i = (M_i) * [sin(θ_i + θM_i), cos(θ_i + θM_i), 0]
        th_m = theta(i) + thetaM(i);
        m_i = M(i) * [sin(th_m); cos(th_m); 0];

        tau_magnetic = m_i(2)*B_ext(1) - m_i(1)*B_ext(2);

        % ---- 3) 스프링 토크 ----
        tau_spring = 0;
        if i > 1
            tau_spring = tau_spring - k_spring(i-1)*(theta(i) - theta(i-1));
        end
        if i < num_links
            tau_spring = tau_spring + k_spring(i)*(theta(i+1) - theta(i));
        end

        % ---- 합력(=0) ----
        eq(i) = tau_magnetic + tau_spring;  
    end
end

%% =========================================================================
%% (D) 링크별 위치 계산 (2D, 순차적으로 이어 붙임)
%% =========================================================================
function positions = compute_link_positions(theta, link_length)
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

