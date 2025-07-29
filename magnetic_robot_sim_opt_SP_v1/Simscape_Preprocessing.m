%% Fixed Condition
% Plate Geometry
% Fixed.PlateLength = 2/1000;
Fixed.PlateLength = 2.35 * 1e-3;
% Fixed.PlateLength = [2.31 2.21 1.89 2.23 2.58 2.58 2.61]
Fixed.PlateWidth = 3/1000;
% Fixed.Thickness = 0.5e-3;
Fixed.Thickness = 0.7e-3;
Fixed.PlateVolume = Fixed.PlateLength*Fixed.PlateWidth*Fixed.Thickness;
Fixed.num_links = 7;
    
% Magnet Properties
Fixed.MagnetHomeConfig = [25/1000 -Fixed.PlateWidth/2];
Fixed.MagnetDirection = [-1 0];
Fixed.Br = 1.0; % [T]
Fixed.Volume = (0.021)^2*(0.01); % [m^3]
Fixed.myu = 4*pi*1e-7;
% Fixed.MagnetMagnitude = Fixed.Br * Fixed.Volume / Fixed.myu;
Fixed.MagnetMagnitude = 4.4942;
Fixed.MangetDipoleMoment = Fixed.MagnetDirection * Fixed.MagnetMagnitude;
    
% Mechanical Properties
Fixed.Kt = 7e-8 * 180/pi; % [Nm/(radian*s)]


%% Simscape Simulation Options
SimParams.gravity_onoff = 0;
SimParams.stiffness = Fixed.Kt/180*pi; % [Nm/(degree*s)]
SimParams.damping = SimParams.stiffness/10;
SimParams.density = 1000;
SimParams.Opacity = 0.5;
SimParams.arrow_Opacity = 1;
SimParams.arrow_scale = 1;
SimParams.force_onoff = 1;
SimParams.torque_onoff = 1;
SimParams.SimulationTime = 1;
SimParams.Prismatic_onoff = 1;
SimParams.ModelName = 'MagneticBasket_Simscape_Optimizer';

open_system(SimParams.ModelName)
set_param(SimParams.ModelName, 'SimMechanicsOpenEditorOnUpdate', 'off');
% set_param('MagneticBasket_Simscape_Optimizer', 'SimulationMode', 'rapid');


%% Input
% Fitness = -4.617
Input.x_init_1 = [0.913521792419827	0.448904782718796	0.432913436158681	0.779105296997416	-0.766830482834701	-0.550594770209900	-0.299881490314832	0.000133355521367087	0.000133546437201994	0.000135197719065930	0.000134779289279190	0.000139826869410360	0.000138777721391194	0.000139646984787483	7.22348734573541e-05];

% Fitness = -4.959
Input.x_init_2 = [1.00916165847214	0.438659822593574	0.494217186344468	0.755306832633380	-0.619764899309504	-0.634045356627496	-0.537836346403770	0.000147086446640267	0.000148170821685519	0.000149972820401047	0.000147571021875325	0.000148656319567183	0.000149789444319739	0.000149362324945304	0.000148527388474000];

% Fitness = -4.417
Input.x_init_test = [deg2rad([55 25 30 45 -45 -35 -20]), repmat(40000, 1, 7) * Fixed.PlateVolume];
Input.x_init_test2 = [deg2rad([0.1 75 65 60 -60 -65 -75]), [10, repmat(40000, 1, 6) 120000] * Fixed.PlateVolume];

%% Configure Functions
addpath("MatlabFunctions")
S.theta = sym('theta', [1,7]); % [radian]
S.DipoleMoment_HomeConfig = sym('DM_HC', [7,2]); % []

S.LinkPoint(1,:) = [S.theta(1) S.theta(1)]*0;
S.theta_acc = 0;
for i = 1:7
    S.theta_acc = S.theta_acc + S.theta(i);
    S.LinkPoint(i+1,:) = S.LinkPoint(i,:) + Fixed.PlateLength * [cos(S.theta_acc), sin(S.theta_acc)];
    S.PlateCenter(i,:) = S.LinkPoint(i,:) + Fixed.PlateLength/2 * [cos(S.theta_acc), sin(S.theta_acc)];
    S.DipoleMoment(i,:) = ([1 0 0;0 1 0] * SymRotation([0 0 1], S.theta_acc) * [S.DipoleMoment_HomeConfig(i,:) 0].').';
end

S.J = jacobian(S.LinkPoint(end,:), S.theta);

Funs.LinkPoint = matlabFunction(S.LinkPoint, 'Vars', {S.theta});
Funs.PlateCenter = matlabFunction(S.PlateCenter, 'Vars', {S.theta});
Funs.DipoleMoment = matlabFunction(S.DipoleMoment, 'Vars', {S.theta, S.DipoleMoment_HomeConfig});



clearvars -except Input Funs SimOptions Fixed SimParams MB