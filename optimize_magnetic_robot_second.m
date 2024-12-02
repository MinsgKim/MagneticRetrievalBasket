function optimize_magnetic_robot_second()
    clc; clear; close;
    % 최적화 파라미터
    num_links = 7; % 고정된 링크 수
    psi_init = 4e04 * ones(1, num_links); % 초기 자화 프로파일 조정
    theta_M_init = rand(1, num_links); % 랜덤 초기화
    r_init = 0.07; % 초기 거리 (단위: 미터)
    link_length_init = 2e-03; % 링크 길이 조정
    disp('Initial theta_M (degrees):');
    disp(rad2deg(theta_M_init))

    % 초기 추측
    x0 = [psi_init, theta_M_init, r_init, link_length_init];

    % 최적화 경계
    lb = [repmat(2e04, 1, num_links), repmat(-pi, 1, num_links), 0.05, 0.001];
    ub = [repmat(6e04, 1, num_links), repmat(pi, 1, num_links), 0.12, 0.005];

    % 고정 위치 제약 조건
    x_fixed_last = 0.0; % 마지막 링크의 고정된 x 위치 조정

    % 외부 자석의 파라미터 설정
    mu0 = 4 * pi * 1e-7; % 진공의 투자율
    external_magnet.Br = 1.22; % 잔류 자화 (단위: 테슬라)
    external_magnet.volume = (0.02)^3; % 자석의 부피 (예: 2cm x 2cm x 2cm 큐브)
    external_magnet.m = external_magnet.Br * external_magnet.volume / mu0; % 자기 모멘트 계산

    % 최적화 결과 저장용 변수 초기화
    x_results = zeros(1, length(x0)); % 반복 횟수를 5로 설정
    cost_values = zeros(1, 1);

    for i = 1:1
        tic;
        % 최적화 옵션
        options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-6, ...
            'ConstraintTolerance', 1e-6, 'MaxFunctionEvaluations', 1e5, ...
            'FiniteDifferenceStepSize', 1e-6, 'OptimalityTolerance', 1e-6, 'Algorithm', 'sqp');

        % 최적화 수행
        x_result = fmincon(@(x) objective_function_mid_position(x, num_links, external_magnet), ...
                            x0, [], [], [], [], lb, ub, ...
                            @(x) nonlcon_position_constraints(x, num_links, x_fixed_last, external_magnet), ...
                            options);

        x_results(i, :) = x_result;
        cost_values(i) = objective_function_mid_position(x_result, num_links, external_magnet);
        disp(['Iteration ', num2str(i), ' completed.']);
        toc;
    end

    % 최소 비용을 갖는 결과 선택
    [~, best_idx] = min(cost_values);
    x_opt = x_results(best_idx, :);

    % 최적화된 파라미터 표시
    disp('최적화된 파라미터:');
    disp(['M: ', num2str(x_opt(1:num_links))]);
    disp(['theta_M (degrees): ', num2str(rad2deg(x_opt(num_links+1:2*num_links)))]);
    disp(['r: ', num2str(x_opt(end-1))]);
    disp(['link_length: ', num2str(x_opt(end))]);

    % 최적화된 파라미터로 로봇 시뮬레이션 및 시각화
    M_opt = x_opt(1:num_links);
    theta_M_opt = x_opt(num_links+1:2*num_links);
    r_opt = x_opt(end-1);
    link_length_opt = x_opt(end);
    M_opt = M_opt * link_length_opt * 0.0033 * 0.0005;

    % 로봇 시뮬레이션
    T_actual_opt = simulate_robot_transform(num_links, M_opt, theta_M_opt, r_opt, link_length_opt, external_magnet);

    % 로봇 시각화
    plot_robot(T_actual_opt);
end

function cost = objective_function_mid_position(x, num_links, external_magnet)
    % 최적화 변수 해제
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    r = x(2*num_links+1);
    link_length = x(end);
    M = M * link_length * 0.0033 * 0.0005;

    % 자기 로봇 시뮬레이션
    T_actual = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet);

    % 중간 링크의 위치 추출
    middle_index = ceil(num_links / 2);
    x_position = T_actual{middle_index}(1:2, 3); % 2D 위치 x 추출
    position_norm = norm(x_position); % 위치의 노름 계산

    % 위치 노름을 최대화하기 위해 음수로 변환하여 반환
    cost = -position_norm;
end

function [c, ceq] = nonlcon_position_constraints(x, num_links, x_fixed_last, external_magnet)
    % 최적화 변수 해제
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    r = x(2*num_links+1);
    link_length = x(end);
    M = M * link_length * 0.0033 * 0.0005;

    % 자기 로봇 시뮬레이션
    T_actual = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet);

    % 제약 조건: 마지막 링크의 x 위치 고정
    last_position = T_actual{end}(1:2, 3); % 마지막 링크의 위치 추출
    ceq = last_position(1) - x_fixed_last; % x 위치는 x_fixed_last와 같아야 함

    % 부등식 제약 조건 없음
    c = [];
end

function T_actual = simulate_robot_transform(num_links, M, theta_M, r, link_length, external_magnet)
    % 파라미터
    k_spring = 0.00092125 * ones(1, num_links-1); % 조인트의 스프링 상수 (PDMS)
    theta_init = zeros(1, num_links); % 초기 링크 각도
    t_span = [0, 10]; % 시간을 늘려서 동역학이 충분히 전개되도록 함
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    % ode45를 사용하여 시뮬레이션
    [~, Y] = ode45(@(t, y) robot_dynamics(t, y, k_spring, M, theta_M, r, link_length, external_magnet), ...
                   t_span, [theta_init, zeros(1, num_links)], options_ode);

    % 변환 행렬 계산
    T_actual = compute_transform_matrices(Y(end, 1:num_links), link_length);
end

function dY = robot_dynamics(~, Y, k_spring, M, theta_M, r, L, external_magnet)
    % 상태 벡터 Y = [theta1, theta2, ..., omega1, omega2, ...]
    num_links = length(M);
    theta = Y(1:num_links);
    omega = Y(num_links+1:end);

    % 가속도 벡터 초기화
    alpha = zeros(num_links, 1);

    % 관성 모멘트 (단위 질량과 길이 L 가정)
    I = (1/12) * L^2;

    % 감쇠 계수
    damping_coefficient = 1e-3;

    % 각 링크의 위치 계산
    positions = compute_link_positions(theta, L);

    % 각 링크에 대한 힘과 토크 계산
    for i = 1:num_links
        % 링크의 위치
        r_vec = positions(:, i) - [-0.00165; num_links*L + r]; % 외부 자석은 robot head 위에 r만큼 위치한다고 가정

        % 외부 자기장 계산
        B_ext = calculate_B_vector(r_vec, external_magnet.m);

        % 링크의 자기 모멘트 벡터
        theta_moment = theta(i) + theta_M(i);
        m_i = M(i) * [sin(theta_moment); cos(theta_moment); 0];

        % 자기 토크 계산 (z축 성분만 고려)
        tau_magnetic = m_i(1) * B_ext(2) - m_i(2) * B_ext(1);

        % 스프링 토크
        tau_spring = 0;
        if i > 1
            tau_spring = tau_spring - k_spring(i-1) * (theta(i) - theta(i-1));
        end
        if i < num_links
            tau_spring = tau_spring + k_spring(i) * (theta(i+1) - theta(i));
        end

        % 감쇠 토크
        tau_damping = -damping_coefficient * omega(i);

        % 모든 토크 합산 및 가속도 계산
        alpha(i) = (tau_magnetic + tau_spring + tau_damping) / I;
    end

    dY = [omega; alpha];
end

function positions = compute_link_positions(theta, link_length)
    % 각 링크의 끝 위치 계산
    num_links = length(theta);
    x = zeros(1, num_links);
    y = zeros(1, num_links);
    current_theta = 0;
    x_current = 0;
    y_current = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        x_current = x_current + link_length * sin(current_theta);
        y_current = y_current + link_length * cos(current_theta);
        x(i) = x_current;
        y(i) = y_current;
    end

    positions = [x; y]; % the end point of each link
end

function T = compute_transform_matrices(theta, link_length)
    % 각 링크에 대한 2D 변환 행렬 계산
    num_links = length(theta);
    T = cell(1, num_links);
    x = 0;
    y = 0;
    current_theta = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        dx = link_length * sin(current_theta);
        dy = link_length * cos(current_theta);

        % 위치 업데이트
        x = x + dx;
        y = y + dy;

        % 변환 행렬 (2D 회전 및 평행 이동)
        T{i} = [cos(current_theta), -sin(current_theta), x; 
                sin(current_theta),  cos(current_theta), y; 
                0,                  0,                  1];
    end
end

function plot_robot(T_actual)
    figure;
    hold on;
    num_links = length(T_actual);
    x_positions = zeros(1, num_links+1);
    y_positions = zeros(1, num_links+1);
    x_positions(1) = 0;
    y_positions(1) = 0;
    for i = 1:num_links
        x_positions(i+1) = T_actual{i}(1,3);
        y_positions(i+1) = T_actual{i}(2,3);
    end
    plot(x_positions, y_positions, '-o', 'LineWidth', 2);
    xlabel('X Position');
    ylabel('Y Position');
    title('Optimized Robot Configuration');
    grid on;
    axis equal;
    hold off;
end

function B = calculate_B_vector(r_vec, m_ext)
    % dipole 모델을 사용하여 외부 자기장 계산
    % r_vec: 위치 벡터 (slave magnet 위치 - 외부 자석 위치)
    % m_ext: 외부 자석의 자기 모멘트 (스칼라 값)

    mu0 = 4 * pi * 1e-7; % 진공의 투자율

    r_norm = norm(r_vec);
    r_hat = r_vec / r_norm;

    % 외부 자석의 자기 모멘트 방향 (y축 방향으로 가정)
    m_vec = [0; m_ext; 0];

    % 위치 벡터를 3D로 확장
    r_vec_3D = [r_vec; 0];

    % 자기장 계산
    B_full = (mu0 / (4 * pi)) * ( (3 * r_vec_3D * (dot(m_vec, r_vec_3D)) / r_norm^5) - (m_vec / r_norm^3) );

    % 2D 평면에서의 자기장 (x, y 성분)
    B = B_full(1:2);
end
