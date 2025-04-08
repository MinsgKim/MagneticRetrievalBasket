% this code is same with optim_bending_and_stability.m, it is written to
% show the y position of each joint

function [Fitness, v] = optim_bending_and_stability_disp(x, SimParams, Options)

arguments

    x (1,:) double
    SimParams (1,1) struct
    Options.DispInfo {mustBeMember(Options.DispInfo, {'yes', 'no'})} = 'no'
    Options.GUI {mustBeMember(Options.GUI, {'on', 'off'})} = 'off'
    Options.InputAngleUnit {mustBeMember(Options.InputAngleUnit, {'Degree', 'Radian'})} = 'Radian'

end

T = tic;

% if angle unit is Degree, change it to Radian
if strcmp(Options.InputAngleUnit, 'Degree')
    x(1:6) = x(1:6) * pi/180;
end

SimParams.Angle = x;

% matlab 작업공간에 'x'라는 변수 넘김
% simulink가 함수 작업공간이 아니라 기본 base 작업공간의 'x'를 참조하기 때문
assignin("base", "Sim", SimParams);

% Simulink simulation
ModelName = 'Optimization_250328';
if strcmp(Options.GUI, 'off')
    set_param(ModelName, 'SimMechanicsOpenEditorOnUpdate', 'off');
else
    set_param(ModelName, 'SimMechanicsOpenEditorOnUpdate', 'on');
end
SimOut = sim(ModelName);
Upper_q = SimOut.Upper_q.Data(end,:);
Lower_q = SimOut.Lower_q.Data(end,:);


    q_acc_upper = 0;
    q_acc_lower = 0;
    LinkPoint_upper(1,:) = [0 0];
    LinkPoint_lower(1,:) = [0 0];
    for i = 1:6
        q_acc_upper = q_acc_upper + Upper_q(i);
        q_acc_lower = q_acc_lower + Lower_q(i);
        LinkPoint_upper(i+1,:) = LinkPoint_upper(i,:) + SimParams.PlateLength * [cos(q_acc_upper), sin(q_acc_upper)];
        LinkPoint_lower(i+1,:) = LinkPoint_lower(i,:) + SimParams.PlateLength * [cos(q_acc_lower), sin(q_acc_lower)];
    end

    F1 = - (max(LinkPoint_upper(:,2))^2+max(LinkPoint_lower(:,2))^2);

    F2 = LinkPoint_upper(end,2)^2 + LinkPoint_lower(end,2)^2;

    w = 0.5;

    Fitness = w*F1+(1-w)*F2;

    v = [LinkPoint_upper(:,2) LinkPoint_lower(:,2)];

end