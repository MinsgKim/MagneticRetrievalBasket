function [Upper_q, Lower_q] = Sim_Optimization(SimParams, Options)
    
    arguments
        SimParams (1,1) struct
        Options.DispInfo {mustBeMember(Options.DispInfo, {'yes', 'no'})} = 'no'
        Options.GUI {mustBeMember(Options.GUI, {'on', 'off'})} = 'off'
    end
    
    T = tic;

    % matlab 작업공간에 'Sim'라는 변수 넘김
    % simulink가 함수 작업공간이 아니라 기본 base 작업공간의 'Sim'를 참조하기 때문
    assignin("base", 'Sim', SimParams);

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

    T = toc;

    % Display Time
    if strcmp(Options.DispInfo, 'yes')
        msg_2 = ['Time : ', num2str(round(T/10000, 2)), '[sec]'];
        disp(msg_2)
    end

end