%% Fixed Condition
% Plate Geometry
Fixed.PlateLength = 2/1000;
Fixed.PlateWidth = 3/1000;
Fixed.PlateThickness = 0.5/1000;
Fixed.PlateVolume = Fixed.PlateLength * Fixed.PlateWidth * Fixed.PlateThickness;
Fixed.num_links = 7;
    
% Magnet Properties
Fixed.AppliedField = 0.015; % [T]
    
% Mechanical Properties
Fixed.Kt = 7e-8 * 180/pi; % [Nm/(radian*s)]


%% Simscape Simulation Options
SimParams.stiffness = Fixed.Kt/180*pi; % [Nm/(degree*s)]
SimParams.damping = SimParams.stiffness/5;
SimParams.density = 2000;

SimParams.SimulationTime = 0.5;
SimParams.ModelName = 'MK_Simscape_Optimizer';

open_system(SimParams.ModelName)
set_param(SimParams.ModelName, 'SimMechanicsOpenEditorOnUpdate', 'off');
% set_param(SimParams.ModelName, 'SimulationMode', 'rapid');


%% Input
% Fitness = -4.617
Input.x_init_1 = [0.913521792419827	0.448904782718796	0.432913436158681	0.779105296997416	-0.766830482834701	-0.550594770209900	-0.299881490314832	0.000133355521367087	0.000133546437201994	0.000135197719065930	0.000134779289279190	0.000139826869410360	0.000138777721391194	0.000139646984787483];

% Fitness = -4.959
Input.x_init_2 = [1.00916165847214	0.438659822593574	0.494217186344468	0.755306832633380	-0.619764899309504	-0.634045356627496	-0.537836346403770	0.000147086446640267	0.000148170821685519	0.000149972820401047	0.000147571021875325	0.000148656319567183	0.000149789444319739	0.000149362324945304];

% Fabricated Model,
Input.x_init_angle= deg2rad([55, 25, 30, 45, -45, -35, -20]);
Input.x_init_magnitude = repmat(40000, 1, Fixed.num_links) * Fixed.PlateVolume;
Input.x_init_test = [Input.x_init_angle, Input.x_init_magnitude];

% random initial parameters
Input.x_init_random_angle = 2*pi*rand(1,Fixed.num_links)-1;
Input.x_init_random_magnitude = (repmat(20000, 1, Fixed.num_links) + 20000 * rand(1,Fixed.num_links)) * Fixed.PlateVolume;
Input.x_init_random = [Input.x_init_random_angle, Input.x_init_random_magnitude];


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