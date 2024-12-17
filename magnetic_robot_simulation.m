classdef magnetic_robot_simulation

    properties
        MagR
        ExMag
    end

    methods
        function obj = magnetic_robot_simulation(Magnetic_Robot, External_Magnet)  % import classes                                 

            obj.MagR = Magnetic_Robot;
            obj.ExMag = External_Magnet;

        end

        function [c, ceq] = nonlcon(obj)

        end

    end

end