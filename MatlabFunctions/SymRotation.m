function Rmatrix = SymRotation(Axis, Angle, Options)
    arguments
        Axis (3,1) double
        Angle (1,1) sym
        Options.AngleUnit {mustBeMember(Options.AngleUnit, {'Degree', 'Radian'})} = 'Radian'
    end

    if strcmp(Options.AngleUnit, 'Degree')
        Angle = Angle/180*pi;
    end

    Axis = Axis/norm(Axis);
    Skew = [0 -Axis(3) Axis(2);Axis(3) 0 -Axis(1);-Axis(2) Axis(1) 0];
    Rmatrix = expm(Skew*Angle);

end