function optimize_magnetic_robot()
    % 최적화 파라미터
    num_links = 7; % 고정된 링크 수
    psi_init = 4e04 * ones(1, num_links); % 초기 자화 프로파일 조정
    theta_M_init = rand(1, num_links); % 제로로 초기화
    B_init = 0.015; % 초기 자기장 증가
    link_length_init = 2e-03; % 링크 길이 조정
    disp(rad2deg(theta_M_init))

    % 초기 추측
    x0 = [psi_init, theta_M_init, B_init, link_length_init];

    % 최적화 경계
    lb = [repmat(2e04, 1, num_links), repmat(-pi, 1, num_links), 0.01, 0.001];
    ub = [repmat(6e04, 1, num_links), repmat(pi, 1, num_links), 0.03, 0.005];

    % 고정 위치 제약 조건
    x_fixed_last = 0.0; % 마지막 링크의 고정된 x 위치 조정

    % 최적화 옵션
    options = optimoptions('fmincon', 'Display', 'iter', 'StepTolerance', 1e-6, ...
        'ConstraintTolerance', 1e-6, 'MaxFunctionEvaluations', 1e4, ...
        'FiniteDifferenceStepSize', 1e-5, 'OptimalityTolerance', 1e-6, 'Algorithm', 'sqp');

    % 최적화 수행
    tic;
    x_opt = fmincon(@(x) objective_function_mid_position(x, num_links), ...
                    x0, [], [], [], [], lb, ub, ...
                    @(x) nonlcon_position_constraints(x, num_links, x_fixed_last), ...
                    options);

    % 최적화된 파라미터 표시
    disp('최적화된 파라미터:');
    disp(['M: ', num2str(x_opt(1:num_links))]);
    disp(['theta_M: ', num2str(rad2deg(x_opt(num_links+1:2*num_links)))]);
    disp(['B: ', num2str(x_opt(end-1))]);
    disp(['link_length: ', num2str(x_opt(end))]);
    toc;
end

function cost = objective_function_mid_position(x, num_links)
    % 최적화 변수 해제
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    B = x(end-1);
    link_length = x(end);
    M = M * link_length * 0.0033 * 0.0005;

    % 외부 자기장을 벡터 형태로 표현
    B_vector = B(1) * [0; -1];

    % 자기 로봇 시뮬레이션
    T_actual = simulate_robot_transform(num_links, M, theta_M, B_vector, link_length);

    % 중간 링크의 변환 행렬에서 위치 추출 (3번째 링크)
    middle_index = ceil(num_links / 2);
    position = T_actual{middle_index}(1:2, 3); % 2D 위치 (x, y) 추출
    position_norm = norm(position); % 위치의 노름 계산

    % 위치 노름을 최대화하기 위해 음수로 변환하여 반환
    cost = -position_norm;
end

function [c, ceq] = nonlcon_position_constraints(x, num_links, x_fixed_last)
    % 최적화 변수 해제
    M = x(1:num_links);
    theta_M = x(num_links+1:2*num_links);
    B = x(end-1);
    link_length = x(end);
    M = M * link_length * 0.0033 * 0.0005;

    % 외부 자기장을 벡터 형태로 표현
    B_vector = B(1) * [0; -1];

    % 자기 로봇 시뮬레이션
    T_actual = simulate_robot_transform(num_links, M, theta_M, B_vector, link_length);

    % 제약 조건: 마지막 링크의 x 위치 고정
    last_position = T_actual{end}(1:2, 3); % 마지막 링크의 위치 추출
    ceq = last_position(1) - x_fixed_last; % x 위치는 x_fixed_last와 같아야 함

    % 부등식 제약 조건 없음
    c = [];
end

function T_actual = simulate_robot_transform(num_links, M, theta_M, B, link_length)
    % 파라미터
    k_spring = 0.00092125 * ones(1, num_links-1);
    theta_init = zeros(1, num_links); % 초기 링크 각도
    t_span = [0, 1]; % 시간 범위 축소
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

    % ode45를 사용하여 시뮬레이션
    [~, Y] = ode45(@(t, y) robot_dynamics(t, y, k_spring, M, theta_M, B, link_length), ...
                   t_span, [theta_init, zeros(1, num_links)], options_ode);

    % 변환 행렬 계산
    T_actual = compute_transform_matrices(Y(end, 1:num_links), link_length);
end

function dY = robot_dynamics(t, Y, k_spring, M, theta_M, B, L)
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

    % 각 링크에 대한 토크 계산
    for i = 1:num_links
        % 자기 토크
        tau_magnetic = M(i) * norm(B) * sin(theta(i) + theta_M(i) - atan2(B(2), B(1)));

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

function T = compute_transform_matrices(theta, link_length)
    % 각 링크에 대한 2D 변환 행렬 계산
    num_links = length(theta);
    T = cell(1, num_links);
    x = 0;
    y = 0;
    current_theta = 0;

    for i = 1:num_links
        current_theta = current_theta + theta(i);
        dx = link_length * cos(current_theta);
        dy = link_length * sin(current_theta);

        % 변환 행렬 (2D 회전 및 평행 이동)
        T{i} = [cos(current_theta), -sin(current_theta), x; 
                sin(current_theta),  cos(current_theta), y; 
                0,                  0,                  1];
        % 위치 업데이트
        x = x + dx;
        y = y + dy;
    end
end