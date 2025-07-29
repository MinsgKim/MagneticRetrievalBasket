%%
clear; clc;
Params = Sim_Param_Setup;

ParamsBus = Simulink.Bus.createObject(Params);

open('tests.slx')


%%
clc;clear;

a = 0.5e-3;
b = 1.5e-3;
c = 2e-3;

D_z = pi^(-1)*((b^2-c^2)/(2*b*c)*log((sqrt(a^2+b^2+c^2)-a)/(sqrt(a^2+b^2+c^2)+a))+...
    (a^2-c^2)/(2*a*c)*log((sqrt(a^2+b^2+c^2)-b)/(sqrt(a^2+b^2+c^2)+b))+b/(2*c)*log((sqrt(a^2+b^2)+a)/(sqrt(a^2+b^2)-a))+...
    a/(2*c)*log((sqrt(a^2+b^2)+b)/(sqrt(a^2+b^2)-b))+c/(2*a)*log((sqrt(b^2+c^2)-b)/(sqrt(b^2+c^2)+b))+...
    c/(2*b)*log((sqrt(a^2+c^2)-a)/(sqrt(a^2+c^2)+a))+2*atan(a*b/c/sqrt(a^2+b^2+c^2))+(a^3+b^3-2*c^3)/(3*a*b*c)+...
    (a^2+b^2-2*c^2)/(3*a*b*c)*sqrt(a^2+b^2+c^2)+c/a/b*(sqrt(a^2+c^2)+sqrt(b^2+c^2))-((a^2+b^2)^(3/2)+(b^2+c^2)^(3/2)+(c^2+a^2)^(3/2))/(3*a*b*c));

a = 2e-3;
b = 0.5e-3;
c = 1.5e-3;

D_z2 = pi^(-1)*((b^2-c^2)/(2*b*c)*log((sqrt(a^2+b^2+c^2)-a)/(sqrt(a^2+b^2+c^2)+a))+...
    (a^2-c^2)/(2*a*c)*log((sqrt(a^2+b^2+c^2)-b)/(sqrt(a^2+b^2+c^2)+b))+b/(2*c)*log((sqrt(a^2+b^2)+a)/(sqrt(a^2+b^2)-a))+...
    a/(2*c)*log((sqrt(a^2+b^2)+b)/(sqrt(a^2+b^2)-b))+c/(2*a)*log((sqrt(b^2+c^2)-b)/(sqrt(b^2+c^2)+b))+...
    c/(2*b)*log((sqrt(a^2+c^2)-a)/(sqrt(a^2+c^2)+a))+2*atan(a*b/c/sqrt(a^2+b^2+c^2))+(a^3+b^3-2*c^3)/(3*a*b*c)+...
    (a^2+b^2-2*c^2)/(3*a*b*c)*sqrt(a^2+b^2+c^2)+c/a/b*(sqrt(a^2+c^2)+sqrt(b^2+c^2))-((a^2+b^2)^(3/2)+(b^2+c^2)^(3/2)+(c^2+a^2)^(3/2))/(3*a*b*c));

a = 1.5e-3;
b = 2e-3;
c = 0.5e-3;

D_z3 = pi^(-1)*((b^2-c^2)/(2*b*c)*log((sqrt(a^2+b^2+c^2)-a)/(sqrt(a^2+b^2+c^2)+a))+...
    (a^2-c^2)/(2*a*c)*log((sqrt(a^2+b^2+c^2)-b)/(sqrt(a^2+b^2+c^2)+b))+b/(2*c)*log((sqrt(a^2+b^2)+a)/(sqrt(a^2+b^2)-a))+...
    a/(2*c)*log((sqrt(a^2+b^2)+b)/(sqrt(a^2+b^2)-b))+c/(2*a)*log((sqrt(b^2+c^2)-b)/(sqrt(b^2+c^2)+b))+...
    c/(2*b)*log((sqrt(a^2+c^2)-a)/(sqrt(a^2+c^2)+a))+2*atan(a*b/c/sqrt(a^2+b^2+c^2))+(a^3+b^3-2*c^3)/(3*a*b*c)+...
    (a^2+b^2-2*c^2)/(3*a*b*c)*sqrt(a^2+b^2+c^2)+c/a/b*(sqrt(a^2+c^2)+sqrt(b^2+c^2))-((a^2+b^2)^(3/2)+(b^2+c^2)^(3/2)+(c^2+a^2)^(3/2))/(3*a*b*c));


clear; clc;

%% ——— 물리 상수 및 시료 형상 ———
clear; clc; close;
mu0 = 4*pi*1e-7;     % [H/m] 자유공간 투자율
Ms  = 1e5;         % [A/m] 포화 자화 (실험값으로 수정)
H   = 2.4/mu0;         % [A/m] 외부장 세기 (실험값으로 수정)
Ku  = 5.6e5;           % [J/m^3] 결정형 anisotropy 상수 (실험값으로 수정)
Ku  = 0;           % [J/m^3] 결정형 anisotropy 상수 (실험값으로 수정)

% 시료 치수 (직육면체 2×1.5×0.5 mm)
L = 2e-3;  W = 1.5e-3;  T = 0.5e-3;
a = L/2;  b = W/2;  c = T/2;

% ——— Aharoni 직육면체 demagnetizing factor 함수 ———
D_prism = @(aa,bb,cc) (1/pi)*((bb^2-cc^2)/(2*bb*cc)*log((sqrt(aa^2+bb^2+cc^2)-aa)/(sqrt(aa^2+bb^2+cc^2)+aa)) ...
    + (aa^2-cc^2)/(2*aa*cc)*log((sqrt(aa^2+bb^2+cc^2)-bb)/(sqrt(aa^2+bb^2+cc^2)+bb)) ...
    + bb/(2*cc)*log((sqrt(aa^2+bb^2)+aa)/(sqrt(aa^2+bb^2)-aa)) ...
    + aa/(2*cc)*log((sqrt(aa^2+bb^2)+bb)/(sqrt(aa^2+bb^2)-bb)) ...
    + cc/(2*aa)*log((sqrt(bb^2+cc^2)-bb)/(sqrt(bb^2+cc^2)+bb)) ...
    + cc/(2*bb)*log((sqrt(aa^2+cc^2)-aa)/(sqrt(aa^2+cc^2)+aa)) ...
    + 2*atan(aa*bb/(cc*sqrt(aa^2+bb^2+cc^2))) ...
    + (aa^3+bb^3-2*cc^3)/(3*aa*bb*cc) ...
    + (aa^2+bb^2-2*cc^2)/(3*aa*bb*cc)*sqrt(aa^2+bb^2+cc^2) ...
    + cc/(aa*bb)*(sqrt(aa^2+cc^2)+sqrt(bb^2+cc^2)) ...
    - (((aa^2+bb^2)^(3/2)+(bb^2+cc^2)^(3/2)+(cc^2+aa^2)^(3/2))/(3*aa*bb*cc)));

% 각 축에 대한 demagnetizing factor
N_x = D_prism(c,b,a);   % easy axis (긴 축) 방향
N_y = D_prism(a,c,b);   % 중간 축
N_z = D_prism(b,a,c);   % 짧은 축

% shape anisotropy 상수 계산 (easy axis를 x축으로 가정)
N_para = N_x;
N_perp = (N_y+N_z)/2;
K_shape = 0.5 * mu0 * Ms^2 * (N_perp - N_para);

% 유효 anisotropy 상수
K_eff = Ku + K_shape;

fprintf('K_shape = %.3e J/m^3,  K_eff = %.3e J/m^3\n', K_shape, K_eff);


% ——— θ-φ 최적화 ———
theta = linspace(0, pi/2, 181);      % 0°~90°
phi_opt = zeros(size(theta));
M_angle = zeros(size(theta));

for i = 1:length(theta)
    th = theta(i);
    % 총 에너지: anisotropy + Zeeman
    E = @(phi) K_eff*sin(th-phi).^2 - mu0*Ms*H*cos(phi);
    phi_opt(i) = fminbnd(E, -pi/2, +pi/2);
    M_angle(i) = theta(i)-phi_opt(i);
end

% ——— 결과 플롯 ———
figure;
plot(theta*180/pi, M_angle*180/pi, 'LineWidth', 2);
grid on;
xlabel('θ (easy‑axis ↔ H) [deg]');
ylabel('θ - φ (easy axis ↔ M) [deg]');
title('φ vs. θ with Effective Anisotropy K_{eff}');
