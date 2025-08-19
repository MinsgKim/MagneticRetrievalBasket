function Params = Sim_Param_Setup

    % Link Geometry (size of each link)
    Params.LinkLength = 1 * 1e-3;
    Params.LinkWidth = 2 * 1e-3;
    Params.LinkThickness = 0.5 * 1e-3;
    Params.LinkVolume = Params.LinkLength * Params.LinkWidth * Params.LinkThickness;
    Params.LinkDensity = 2000;
    
    % Link Magnetic Moment (set magnetic magnitude & direction of each link)
    Params.Link_Magnetization_Magnitude = 40e3*ones(1,6);
    Params.Link_Magnetization_Angle = [25 30 45 -45 -35 -25]; % degree
    Params.Link_Moment = Params.Link_Magnetization_Magnitude.*Params.LinkVolume;
    
    % Source Magnet Properties (external magnetic source)
    Params.Source_MagnetMagnitude = 0;
    % Params.Source_MagnetMagnitude = 4.4942;
    Params.Source_Angle = 0;
    Params.Source_Distance_x = 60e-3;
    Params.Source_Distance_y = 15e-3;
    
    % Mechanical Properties (spring constant, damping coefficient, vacuum permeability)
    Params.Stiffness = 0.9e-7; % [Nm/(radian)]
    Params.Damping = 1.5e-8; % [Nm/(radian/s)]
    Params.mu = 4*pi*1e-7;

    % etc (divide a single link into multiple minute links for estimating
    % a better moment)
    Params.N = 16; % N x N for finite element of adjacent link
    
end