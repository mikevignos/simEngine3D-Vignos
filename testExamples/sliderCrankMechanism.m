%% sliderCrankMechanism.m
% Taken from ch. 12.2 section Edward J. Haug: Computer Aided Kinematics and Dynamics of 
% Mechanical Systems (Allyn and Bacon, 1989)
% Used to test the combinaton of a revolute joint and translational joint
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 1;
length1 = 0;
isGround1 = 1;
JMatrix1 = eye(3,3);
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1);

% Add body2. The is the crank
length2 = 0.08; % meters
mass2 = 0.12; % kg

% Define inertial properties of the crank
Jxx = 0.0001;
Jyy = 0.00001;
Jzz = 0.0001;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2);

% Add body3 to the system. This is the connector between the slider and the
% crank.
length3 = 0.3; % meters
mass3 = 0.5; % kg

% Define inertial properties of the connecting rod
Jxx = 0.004;
Jyy = 0.0004;
Jzz = 0.004;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3);

% Add the slider (body4) to the mechanism.
length4 = 0.05; % The length of the slider does not matter so this is an arbitrary length.
mass4 = 2.0;

% Define inertial properties of the slider
Jxx = 0.0001;
Jyy = 0.0001;
Jzz = 0.0001;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'block', isGround4, mass4, length4, JMatrix4);





%%%%%%%%%%%%%%%%%%%% For this example, I need to implement the ability to
%%%%%%%%%%%%%%%%%%%% prescribe a constant angular velocity.


%% Set initial conditions of each body
% The initial positions and orientations are provided from section 10.2.2.
% Initial position
% r1Initial = zeros(3,1); % Ground
% r2Initial = [0.00002; 0.09982; 0.12005]; % Crank
% r3Initial = [0.09993; 0.05183; 0.09998]; % Connector
% r4Initial = [0.19959; 0.00057; -0.00008]; % Slider.
% rInitial = [r1Initial; r2Initial; r3Initial; r4Initial];
% 
% % Initial orientation
% % Ground
% p1 = [1.0 0.0 0.0 0.0]'; 
% 
% % Crank
% p2 = [0.72090 0.69306 0.00004 -0.00004]';
% 
% % Connector (body3)
% p3 = [0.88723 -0.21202 0.39833 -0.09569]';
% 
% % Slider
% p4 =  [1.0 0.0 -0.00012 -0.00025]';
% 
% pInitial = [p1; p2; p3; p4];

r1Initial = zeros(3,1); % Ground
r2Initial = [0.0; 0.1; 0.12]; % Crank
r3Initial = [0.1; 0.05; 0.1]; % Connector
r4Initial = [0.2; 0.0; 0.0]; % Slider.
rInitial = [r1Initial; r2Initial; r3Initial; r4Initial];

% Initial orientation
% Ground
p1 = [1.0 0.0 0.0 0.0]'; 

% Crank
p2 = [0.7042 0.71 0.0 0.0]';

% Connector (body3)
p3 = [0.8865 -0.21 0.4 -0.1]';

% Slider
p4 =  [1.0 0.0 0.0 0.0]';

pInitial = [p1; p2; p3; p4];

% Initial velocities
% rDotInitial = zeros(6,1);
% pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, [], [], pInitial, [], [], t);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])
% saveas(gcf,'A8P1_MechanismInitialPosition.png');

%% Define revolute joint between crank and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0.1 0.12]';
a1.pointOnBody2 = [0 0 0]';
a1.vector1OnBody1 = [0 1 0]';
a1.vector2OnBody1 = [0 0 1]';
a1.vectorOnBody2 = [1 0 0]';
a1.constraintName = 'Revolute Joint b/w Ground and Crank';

sys.addJoint('revolute',a1);

%% Define spherical joint between connector and crank
a2.body1 = 2;
a2.body2 = 3;
a2.pointOnBody1 = [0.0 0.08 0.0]';
a2.pointOnBody2 = [-0.15 0.0 0.0]';  
a2.constraintName = 'Spherical joint b/w connector and crank';

sys.addJoint('spherical',a2);

%% Define revolute-cylindrical composite joint between connector and slider
a3.body1 = 3;
a3.body2 = 4;
a3.pointOnBody1 = [0.15 0 0]';
a3.pointOnBody2 = [1 0 0]';
a3.vectorOnBody1 = [0 1 0]';
a3.vector1OnBody2 = [1 0 0]';
a3.vector2OnBody2 = [0 1 0]';
a3.vector3OnBody2 = [0 0 1]';
a3.constraintName = 'Rev-cylindrical joint b/w connector and slider';

sys.addJoint('revolute-cylindrical',a3);

%% Define translation joint between slider and ground
a4.body1 = 1;
a4.body2 = 4;
a4.pointOnBody1 = [1 0 0]';
a4.pointOnBody2 = [0 0 0]';
a4.vector1OnBody1 = [0 1 0]';
a4.vector2OnBody1 = [0 0 1]';
a4.vector1OnBody2 = [1 0 0]';
a4.vector2OnBody2 = [0 0 1]';
a4.constraintName = 'Translation joint b/w slider and ground';

sys.addJoint('translational',a4);

%% Define distance constraint between connector and slider
necessaryAttributes = [{'bodyI'} {'bodyJ'} {'sBarIP'} {'sBarJQ'} {'ft'} {'ftDot'} {'ftDDot'}];
a5.bodyI = 3;
a5.bodyJ = 4;
a5.sBarIP = [0.15 0 0]';
a5.sBarJQ = [0 0 0]';
a5.ft = @(t)0;
a5.ftDot =  @(t)0;
a5.ftDDot =  @(t)0;
a5.constraintName = 'Distance constraint b/w connector and slider';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'d',a5);

%% Add driving constraint to model
% DP1 constraint between -Z and y'
a6.bodyJ = 1;
a6.bodyI = 2;
a6.aBarJ = [0 1 0]';
a6.aBarI = [0 0 -1]';
a6.ft = @(t)cos(-2*pi*t);
a6.ftDot = @(t)(2*pi*sin(-2*pi*t));
a6.ftDDot = @(t)(4*pi^2*cos(-2*pi*t));
a6.constraintName = 'DP1 driving constraint';
isKinematic = 0;
sys.addBasicConstraint(isKinematic,'dp1',a6);

%% Perform kinematics analysis
if 1
    timeStart = 0;
    timeEnd = 1;
    timeStep = 10^-3;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
    sys.inverseDynamicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
%     sys.kinematicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
    analysisTime = toc;
    save('sliderCrankMechanism.mat','sys');
else
    load('sliderCrankMechanism.mat')
end

% disp(['Inverse Dynamics Analysis for sliderCrankMechanism took ' num2str(analysisTime) ' seconds.'])

%% Plot the position of the slider vs time
sliderPosition = sys.myBodies{4}.myRTotal;
time = sys.myBodies{4}.myTimeTotal;

figure
hold on
plot(time,sliderPosition(1,:))
plot(time,sliderPosition(2,:))
plot(time,sliderPosition(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Position (m)')
title('Slider')
hold off

%% Plot the position of point B versus time
crankOrientation = sys.myBodies{2}.myPTotal;
crankPosition = sys.myBodies{2}.myRTotal;
pointBPosition = zeros(3,length(time));
sB = [0 0.08 0]';
for iT = 1:length(time)
    A = simEngine3DUtilities.p2A(crankOrientation(:,iT));
    pointBPosition(:,iT) = crankPosition(:,iT) + A*sB;
end

figure
hold on
plot(time,pointBPosition(1,:));
plot(time,pointBPosition(2,:));
plot(time,pointBPosition(3,:));
legend('x','y','z')
hold off

%% Display torque at the revolute joint 
% Extract torque for body 2 due to all constraints and time
torque = sys.myBodies{2}.myConstraintTorquesOmegaTotal;
time = sys.myBodies{2}.myTimeTotal;

% Compute theta over time.
theta = 2*pi*time + pi/2;

% Extract torque due to DP1 driving constraint.
DP1const = 18;
torqueDriving = torque((3*DP1const-2):3*DP1const,:);

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
plot(time,torqueDriving(3,:));
xlabel('Time (sec)')
ylabel('Torque (N*m)')
% axis([0 1 -1 1]);
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ')
title('Torque Due to DP1 Driving Constraint in Pendulum Reference Frame')
saveas(gcf,'A8P1_TorqueVsTime.png')

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
[ax, h1, h2] = plotyy(time,torqueDriving(3,:),time,theta);
ax(1).YLim = [-250 250];
xlabel('Time (sec)')
ylabel(ax(1),'Torque (N*m)')
ylabel(ax(2),'Theta (rad)')
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ','Theta')
saveas(gcf,'A8P1_TorqueAndThetaVsTime.png')