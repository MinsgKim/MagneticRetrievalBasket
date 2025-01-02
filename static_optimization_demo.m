function static_optimization_demo()
    %% 1. 기본 파라미터 설정
    clc; clear; close all;

    % 로봇 링크 개수
    num_links = 7;
    % 링크 길이 (m)
    link_length = 2e-3;        % 2 mm
    % 자화(psi) 단위 변환용 단면적 (m^2)
    cross_section_area = 0.0033 * 0.0005;  % (3.3 mm x 0.5 mm)
    % 스프링 계수(링크 사이)
    %  -> 매우 작게 해서 잘 휘도록 설정 (예: 1e-5 N·m/rad)
    k_spring = 7e-8 * ones(1, num_links-1);

    % 외부 자석과의 거리(대략) - Dipole 모델에서 참조
    r_ext = 0.05;  % 40 mm (예시)
    % 외부 자석 객체 (이미 클래스가 있다고 가정)
    EM = External_Magnet();

    %--------------------------------------------
    % 최적화 변수: x = [psi_1..psi_n,  thetaM_1..thetaM_n]
    %--------------------------------------------
    % 자화 세기 범위 (A/m)
    lb_psi = 2e4;  ub_psi = 6e4;
    % 자화 방향 범위 (rad)
    lb_thetaM = -pi;  ub_thetaM = pi;

    % 초기값
%     rng(0);  % 재현성
    psi_init    = 4e4 * (rand(1,num_links) + 0.5);  % [2e4~6e4 사이 랜덤]
    thetaM_init = (rand(1,num_links)*2*pi - pi);    % [-pi~+pi 사이 랜덤]
    x0 = [psi_init, thetaM_init];

    % 하한/상한
    lb = [lb_psi*ones(1,num_links),    lb_thetaM*ones(1,num_links)];
    ub = [ub_psi*ones(1,num_links),    ub_thetaM*ones(1,num_links)];

    %--------------------------------------------
    % fmincon 옵션 + GlobalSearch 예시
    %--------------------------------------------
%     opts = optimoptions('fmincon',...
%         'Display','iter',...
%         'MaxFunctionEvaluations',1e5,...
%         'MaxIterations',1000,...
%         'Algorithm','interior-point',...
%         "EnableFeasibilityMode",true,...
%         "SubproblemAlgorithm","cg");

    opts = optimoptions('fmincon',...
        'Display','iter',...
        'MaxFunctionEvaluations',1e5,...
        'MaxIterations',1000,...
        'Algorithm','sqp');

    problem = createOptimProblem('fmincon',...
        'x0', x0,...
        'lb', lb, 'ub', ub,...
        'objective', @(x) objective_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
        'nonlcon',  @(x) constraint_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM),...
        'options', opts);

    gs = GlobalSearch;
    [x_opt, fval_opt] = run(gs, problem);

    fprintf('==== Optimization Finished ====\n');
    fprintf('Best Cost: %.4f\n', fval_opt);
    disp('Optimized Variables x_opt = [psi_1..psi_n, thetaM_1..thetaM_n]:');
    disp(x_opt);
    

    %% 2. 최적해 시뮬레이션 (정적 해석)
    psi_opt    = x_opt(1:num_links);
    thetaM_opt = x_opt(num_links+1 : 2*num_links);

    % 최적 자화 -> m (A·m^2)로 변환
    M_opt = psi_opt .* (link_length * cross_section_area);

    % 정적 해석으로 링크 각도(theta)를 구함
    theta_final = solve_static_equilibrium(num_links, M_opt, thetaM_opt, r_ext, link_length, EM, k_spring);

    % 링크 위치 계산
    positions_opt = compute_link_positions(theta_final, link_length);

    %% 3. 결과 출력 및 간단 플롯
    figure; hold on; axis equal;

    for i = 1:num_links
        r_vec = positions_opt(:,i) - [0; r_ext];  % 혹은 자석 위치
        B_local = EM.Cal_B(r_vec);
        fprintf('Link %d:  r=%.4e, |B|=%.4e\n', i, norm(r_vec), norm(B_local));
    end


    plot([0, positions_opt(1,1)], [0, positions_opt(2,1)], 'ko-','LineWidth',2);  % 첫 링크 연결

    for i = 1:num_links-1
        plot(positions_opt(1,i:i+1), positions_opt(2,i:i+1), 'bo-','LineWidth',2);
    end
    plot(positions_opt(1,:), positions_opt(2,:), 'ro','MarkerSize',8,'LineWidth',2);

    title('Final Robot Configuration (Static Equilibrium)');
    xlabel('X (m)');
    ylabel('Y (m)');
    grid on;

    % 시각적으로 링크 자화 방향 표시(간단 예)
    current_th = 0;
    for i = 1:num_links
        current_th = current_th + theta_final(i);
        % 자화 각도
        th_m = theta_final(i) + thetaM_opt(i);
        % 가운데 지점
        if i == 1
            x_center = positions_opt(1,i)/2;
            y_center = positions_opt(2,i)/2;
        else
            x_center = (positions_opt(1,i)+positions_opt(1,i-1))/2;
            y_center = (positions_opt(2,i)+positions_opt(2,i-1))/2;
        end
        quiver(x_center, y_center, ...
            0.5*link_length*sin(th_m), ...
            0.5*link_length*cos(th_m), ...
            'Color',[1,0,0],'LineWidth',1.5,'MaxHeadSize',2);
    end

    hold off;

end

%% =========================================================================
%% (A) 목적함수: 중간 링크들의 x좌표 합을 "최대한" 크게 -> 음수 부호로 최소화
%% =========================================================================
function cost = objective_static(x, num_links, link_length, cross_section_area, r_ext, k_spring, EM)
    % x = [psi_1..psi_n, thetaM_1..thetaM_n]
    psi = x(1:num_links);
    thetaM = x(num_links+1 : 2*num_links);

    % A/m -> 실제 자기모멘트 크기(A·m^2)
    M = psi .* (link_length * cross_section_area);

    % 정적 해석
    theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring);

    % 링크 위치 계산
    pos = compute_link_positions(theta_eq, link_length);

    % 중간 링크(2~(n-1))의 x좌표 합
    if num_links > 2
        x_mid = pos(1, 2 : end-1);
    else
        % 링크가 2개 이하라면 중간 링크가 없음
        x_mid = 0;
    end

    % 최대화 -> -sum()
    cost = -sum(x_mid);
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

    x_con1 = x_end - 0.001;   % <=0  -> x_end <= +1 mm
    x_con2 = -x_end - 0.001;  % <=0  -> x_end >= -1 mm
%     y_con  = 0.003 - y_end;   % <=0  -> y_end >= 3 mm

    % 등가 제약(모두 0이어야 함)
%     ceq = [ceq_first_x; ceq_first_y];
%     ceq = [ceq_first_y];

    % 부등호 제약(모두 <= 0)
%     c = [c_link_dist(:); x_con1; x_con2];
end

%% =========================================================================
%% (C) 정적 해석(각 링크에 토크합=0) -> fsolve
%% =========================================================================
function theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring)
    % 초기값 (약간 랜덤 or 0)
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

