function main_optimization_fromFlat()
    clear; clc; close all;

    %% 1) 파라미터
    num_links = 7;
    link_length = 2e-3;
    cross_section_area = 0.0033*0.0005;
    r = 30e-3;  % 외부자석 x위치 (약 30mm)
    k_spring = 1e-05 * ones(1,num_links);

    % 초기 추정값 설정:
    %  (a) psi(자화 세기), theta_M(자화 각도)는 임의/랜덤
    %  (b) theta_link = 0(모두 x축에 누워 있다)
    psi_init = 5000 * (7*rand(1,num_links)+1);
    thetaM_init = pi/4 * (2*rand(1,num_links)-1);
    thetaLink_init = zeros(1,num_links);  % "누워있는" 상태
    x_init = [psi_init, thetaM_init, thetaLink_init];

    % 상하한 설정
    lb = [repmat(5e03, 1, num_links), ...
          repmat(-pi,   1, num_links), ...
          repmat(-pi/2,   1, num_links)];
    ub = [repmat(6e04, 1, num_links), ...
          repmat(pi,    1, num_links), ...
          repmat(pi/2,    1, num_links)];

    %% 2) 객체들 생성
    cf = cost_function();
    EM2 = External_Magnet2();
    RS = RobotState();

    % fmincon 옵션
    options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',1e5,...
        'Algorithm','interior-point',"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg");

    problem = createOptimProblem('fmincon',...
        'x0', x_init, ...
        'lb', lb, 'ub', ub, ...
        'objective', @(x) cf.moment_equilibrium(x, num_links, link_length, cross_section_area, r, k_spring, EM2),...
        'nonlcon', @(x) cf.nonlcon(x, num_links, link_length, cross_section_area, r, k_spring, EM2),...
        'options', options);

    gs = GlobalSearch;  % or MultiStart
    [x_opt, fval_opt] = run(gs, problem);

    fprintf("==== Optimization Finished ====\n");
    fprintf("Best Cost (fval_opt): %.4f\n", fval_opt);

    %% 3) 결과 파싱
    psi_opt = x_opt(1:num_links);
    thetaM_opt = x_opt(num_links+1:2*num_links);
    thetaL_opt = x_opt(2*num_links+1:end);

    disp("Optimized psi:");
    disp(psi_opt);
    disp("Optimized thetaM (deg):");
    disp(rad2deg(thetaM_opt));
    disp("Optimized thetaLink (deg):");
    disp(rad2deg(thetaL_opt));

    %% 4) 토크, 최종 형상 확인
    tau_opt = cf.get_tau(x_opt, num_links, link_length, cross_section_area, r, k_spring, EM2);
    disp("Tau (mag & spring):");
    disp(tau_opt);

    % 시각화
    figure; hold on;
    title('Final Robot Configuration');
    joint_pos = RS.set_joint_pos(num_links, link_length, thetaL_opt);
    plot(joint_pos(1,:), joint_pos(2,:), '-o','LineWidth',2);
    axis equal; grid on; xlabel('X'); ylabel('Y');

    link_center = RS.set_link_center(num_links, link_length, thetaL_opt);
    theta_sum = 0;
    for i=1:num_links
        theta_sum = theta_sum + thetaL_opt(i);
        tM_abs = theta_sum + thetaM_opt(i);
        xMag = cos(tM_abs);
        yMag = sin(tM_abs);
        quiver(link_center(1,i), link_center(2,i), xMag, yMag, 0.002, 'r','LineWidth',1.5);
    end

    disp("Configuration plotting done.");
end
