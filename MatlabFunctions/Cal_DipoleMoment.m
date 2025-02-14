function [DipoleMoment] = Cal_DipoleMoment(Angle_Magnitude, Options)

    arguments
        Angle_Magnitude (1,:) double
        Options.AngleUnit {mustBeMember(Options.AngleUnit, {'Degree', 'Radian'})} = 'Degree'
    end
    
    Angle = Angle_Magnitude(1:7);
    Magnitude = Angle_Magnitude(8:end);

    if strcmp(Options.AngleUnit, 'Degree')
        theta_r = Angle * pi/180;
    else
        theta_r = Angle;
    end

    for i = 1:length(Angle)
        R = axang2rotm([0 0 1 theta_r(i)]);
        DipoleMoment(i,:) = ([1 0 0;0 1 0]*R*[1 0 0].').' * Magnitude(i);
    end

end