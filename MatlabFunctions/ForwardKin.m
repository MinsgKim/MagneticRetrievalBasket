function [LinkPoint, PlateCenter, DipoleMoment] = ForwardKin(q, Angle_Magnitude, Fixed, Options)

    arguments
        q (7,1) double
        Angle_Magnitude (1,:) double
        Fixed (1,1) struct
        Options.Plot {mustBeMember(Options.Plot, {'yes', 'no'})} = 'no'
    end

    DipoleMoment_HomeConfig = Cal_DipoleMoment(Angle_Magnitude);
    q_acc = 0;
    LinkPoint(1,:) = [0 0];
    for i = 1:7
        q_acc = q_acc + q(i);
        LinkPoint(i+1,:) = LinkPoint(i,:) + Fixed.PlateLength * [cos(q_acc), sin(q_acc)];
        PlateCenter(i,:) = LinkPoint(i,:) + Fixed.PlateLength/2 * [cos(q_acc), sin(q_acc)];
        DipoleMoment(i,:) = ([1 0 0;0 1 0] * axang2rotm([0 0 1 q_acc]) * [DipoleMoment_HomeConfig(i,:) 0]')';
    end

end