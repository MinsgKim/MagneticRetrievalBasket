function test_equilibrium()
    % 간단한 테스트 스크립트
    % 1) 파라미터/객체 설정
    % 2) solve_static_equilibrium 호출
    % 3) 결과(각도, 좌표) 확인
    
    clc; clear; close all;

    % ----- 1) 로봇 / 자석 파라미터 설정 -----
    num_links         = 7;
    link_length       = 2e-3;                % 2 mm
    cross_section_area= 0.0033 * 0.0005;     % (3.3 mm x 0.5 mm)
    r_ext            = 0.045;               % 25 mm (자석 위치)
    k_spring         = 7e-8 * ones(1,num_links);  % 매우 작은 스프링 계수
    
    % 외부 자석 객체 (같은 폴더 내 External_Magnet.m 필요)
    EM = External_Magnet(); 
    
    % ----- 2) 임의의 psi, thetaM 값 설정 -----
    psi_test    = 4.0e4 * ones(1, num_links);  % 자화 세기
%     thetaM_test = 2*pi*rand(1, num_links) - pi;         % 자화 방향 = 0
    thetaM_test = [pi, pi, pi, pi, -pi, -pi, -pi];
    disp(thetaM_test)
    M_test      = psi_test .* (link_length * cross_section_area);
    
    % ----- 3) solve_static_equilibrium 호출 -----
    theta_eq_test = solve_static_equilibrium(num_links, ...
                        M_test, thetaM_test, r_ext, link_length, EM, k_spring);
    
    % ----- 4) 링크 위치 계산 및 결과 출력 -----
    positions_test = compute_link_positions(theta_eq_test, link_length);

    for i = 1:num_links
        r_vec = positions_test(:,i) - [0; r_ext];  % 혹은 자석 위치
        B_local = EM.Cal_B(r_vec);
        fprintf('Link %d:  r=%.4e, |B|=%.4e\n', i, norm(r_vec), norm(B_local));
    end
    
    disp('===== 결과 확인 =====');
    disp('최종 theta_eq_test ='); 
    disp(theta_eq_test);
    disp('링크 끝단 좌표 positions_test ='); 
    disp(positions_test);
end


%% ========================================================================
%%  (A) 정적 해석: 자화 토크 + 스프링 토크 = 0 인 각도 벡터 theta 찾기
%% ========================================================================
function theta_eq = solve_static_equilibrium(num_links, M, thetaM, r_ext, link_length, EM, k_spring)
    % fsolve를 이용해, 각 링크별 토크합 = 0인 theta를 찾는다.

    % 초기 각도 (아주 작은 랜덤값)
    theta_init = 1e-2 * randn(1, num_links);

    opts = optimoptions('fsolve', ...
        'Display','none', ...
        'MaxIterations',1000, ...
        'MaxFunctionEvaluations',1e5, ...
        'FunctionTolerance',1e-12);

    theta_eq = fsolve(@(theta) equilibrium_equations(theta, num_links, M, thetaM, ...
                                  r_ext, link_length, EM, k_spring), ...
                      theta_init, opts);
end


%% ========================================================================
%%  (B) 각 링크별 토크합(=0) 방정식
%% ========================================================================
function eq = equilibrium_equations(theta, num_links, M, thetaM, r_ext, link_length, EM, k_spring)
    % theta : 1 x num_links
    % 각 링크 i에 대해 토크 = 자화토크 + 스프링토크 = 0
    
    eq = zeros(1, num_links);

    % 링크 i 끝단의 (x,y) 좌표
    positions = compute_link_positions(theta, link_length);

    for i = 1:num_links
        % -- 1) 외부 자기장 계산 --
        %    자석이 (0, r_ext)에 있다고 가정
        r_vec = positions(:, i) - [0; r_ext];
        B_ext = EM.Cal_B(r_vec);

        % -- 2) 자화토크 (2D z성분) --
        th_m = theta(i) + thetaM(i);
        m_i = M(i)*[sin(th_m); cos(th_m); 0];  % 3D 벡터 (z=0)
        tau_magnetic = m_i(2)*B_ext(1) - m_i(1)*B_ext(2);

        % -- 3) 스프링 토크 --
        tau_spring = 0;
        if i > 1
            tau_spring = tau_spring - k_spring(i-1)*(theta(i) - theta(i-1));
        end
        if i < num_links
            tau_spring = tau_spring + k_spring(i)*(theta(i+1) - theta(i));
        end

        % -- 4) 토크 합계 --
        eq(i) = tau_magnetic + tau_spring;
    end
end


%% ========================================================================
%%  (C) 링크별 좌표 계산
%% ========================================================================
function positions = compute_link_positions(theta, link_length)
    % theta(i)를 순차적으로 누적해 (x,y)를 구한다
    n = length(theta);
    positions = zeros(2, n);

    x = 0; 
    y = 0;
    current_angle = 0;
    for i = 1:n
        current_angle = current_angle + theta(i);
        x = x + link_length * sin(current_angle);
        y = y + link_length * cos(current_angle);
        positions(:, i) = [x; y];
    end
end
