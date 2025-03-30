function MB = Call_MagneticBasket_Params

    % Plate Geometry
    MB.PlateLength = 2.35 * 1e-3;
    MB.PlateWidth = 3.0 * 1e-3;
    MB.PlateThickness = 0.7 * 1e-3;
    MB.PlateVolume = MB.PlateLength * MB.PlateWidth * MB.PlateThickness;
    MB.PlateDensity = 2000;
        
    % Plate Magnet Properties
    % MB.Plate_MagnetMagnitude = [40000 40000 40000 40000 40000 40000] * MB.PlateVolume;
    % magnetization-based magnitude
    MB.Plate_MagnetMagnitude = ones(1,6) * 80000 * MB.PlateVolume;
    % dipole-moment-based magnitude
    % MB.Plate_MagnetMagnitude = ones(1,6) * 15 *1e-4;


    % Source Magnet Properties
    MB.BField_Magnitude = 20;
    MB.Source_MagnetMagnitude = 4.4942; 

    % Mechanical Properties
    MB.Stiffness = 7e-8; % [Nm/(radian)]
    MB.Damping = 5e-9*2; % [Nm/(radian/s)]
    MB.Stiffness_upper = 7e-8*2;
    MB.Damping_upper = 5e-9*2;
    MB.Transition_upper = 0.01;

    % Simulink Parameters
    MB.Opt_EndTime = 0.15; % [sec]
    MB.Cal_EndTime = 0.5; % [sec]
    
end