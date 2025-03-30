function cost = calibrate_Kt(x, q_measured, SimParams, Options)

arguments

    x (1,:) double  % [dipole moment of each link (1x6), Kt]
    q_measured (1,:) double
    SimParams (1,1) struct
    Options.DispInfo {mustBeMember(Options.DispInfo, {'yes', 'no'})} = 'no'
    Options.GUI {mustBeMember(Options.GUI, {'on', 'off'})} = 'off'

end

% T = tic;

SimParams.Stiffness = x(end-1);
SimParams.Stiffness_upper = x(end);
SimParams.Plate_MagnetMagnitude = x(1:6);

% matlab 작업공간에 'Sim'라는 변수 넘김
% simulink가 함수 작업공간이 아니라 기본 base 작업공간의 'Sim'를 참조하기 때문
assignin("base", 'Sim', SimParams);

% Simulink simulation
ModelName = 'Calibration_250328';
if strcmp(Options.GUI, 'off')
    set_param(ModelName, 'SimMechanicsOpenEditorOnUpdate', 'off');
else
    set_param(ModelName, 'SimMechanicsOpenEditorOnUpdate', 'on');
end
SimOut = sim(ModelName);
q = SimOut.q.Data(end,:);

% T = toc;

cost = norm(abs(-q - q_measured));

end