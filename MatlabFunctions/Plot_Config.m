function Plot_Config(Angle_Magnitude, q, Fixed, Options)

    arguments
        Angle_Magnitude (1,14) double
        q (1,7) double
        Fixed (1,1) struct
        Options.ArrowScale (1,1) double = 1
        Options.LineStyle {mustBeMember(Options.LineStyle, {'-', '--'})} = '-'
    end
    
    [LinkPoint, PlateCenter, DipoleMoment] = ForwardKin(q, Angle_Magnitude, Fixed);
    
    plot(LinkPoint(:,1), LinkPoint(:,2), 'LineStyle', Options.LineStyle, 'Marker', '*', 'LineWidth', 1.5);
    hold on
    axis equal
    grid on
    scale = Options.ArrowScale * 25;
    quiver(PlateCenter(:,1), PlateCenter(:,2), DipoleMoment(:,1)*scale, DipoleMoment(:,2)*scale, 'LineWidth', 1.5, 'AutoScale', 'off')
    xlim([-2 16]/1000)
    ylim([-3 7]/1000)
    set(gcf, 'Color', 'w');

end