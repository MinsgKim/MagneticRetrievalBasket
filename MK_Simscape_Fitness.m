function Fitness = MK_Simscape_Fitness(x, SimParams, Fixed, Options)
    
    arguments
        x (1,:) double
        SimParams (1,1) struct
        Fixed (1,1) struct
        Options.DispInfo {mustBeMember(Options.DispInfo, {'yes', 'no'})} = 'no'
        Options.InputAngleUnit {mustBeMember(Options.InputAngleUnit, {'Degree', 'Radian'})} = 'Radian'
        Options.GUI {mustBeMember(Options.GUI, {'on', 'off'})} = 'off'
    end
    
    T = tic;

    % if angle unit is Degree, change it to Radian
    if strcmp(Options.InputAngleUnit, 'Degree')
        x(1:7) = x(1:7) * pi/180;
    end

    % matlab 작업공간에 'x'라는 변수 넘김
    % simulink가 함수 작업공간이 아니라 기본 base 작업공간의 'x'를 참조하기 때문
    assignin("base", 'x', x);

    % Simulink simulation
    if strcmp(Options.GUI, 'off')
        set_param(SimParams.ModelName, 'SimMechanicsOpenEditorOnUpdate', 'off');
    else
        set_param(SimParams.ModelName, 'SimMechanicsOpenEditorOnUpdate', 'on');
    end
    SimOut = sim(SimParams.ModelName);
    Sim_q = SimOut.Config.Data(end,:)*pi/180;

    [LinkPoint, ~, ~] = ForwardKin(Sim_q, x, Fixed);
    Fitness = -max(abs(LinkPoint*[0;1]));
    T = toc;

    if strcmp(Options.DispInfo, 'yes')
        msg_1 = ['Time : ', num2str(round(T/1000,2)), '[sec]'];
        msg_2 = ['Fitness : ', num2str(round(Fitness*1000, 3)), '[mm]'];
        disp(' ')
        disp(msg_1)
        disp(msg_2)
        disp('---------------------')
    end

end