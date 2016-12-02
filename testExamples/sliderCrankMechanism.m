%% sliderCrankMechanism.m
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
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
% Initial position
r1Initial = zeros(3,1); % Ground
r2Initial = [0; 0.1; 0.12]; % Crank
r3Initial = []; % Connector
r4Inital = []; % Slider.
rInitial = [r1Initial; r2Initial];

% Initial orientation
p1 = [1 0 0 0]'; 
s2 = sqrt(2)/2;
A2 = [0 0 1;
    s2 s2 0;
    -s2 s2 0];
p2 = simEngine3DUtilities.A2p(A2);
pInitial = [p1; p2];

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Plot starting configuration of system
% sys.plot(1);
% view([90 0])
% saveas(gcf,'A8P1_MechanismInitialPosition.png');

%% Define revolute joint
necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} {'pointOnBody2'} {'vector1OnBody1'} {'vector2OnBody1'} {'vectorOnBody2'}];
a.body1 = 1;
a.body2 = 2;
a.pointOnBody1 = [0 0 0]';
a.pointOnBody2 = [-2 0 0]';
a.constraintName = 'Spherical Joint';

sys.addJoint('spherical',a);

%% Add driving constraint to model
% % DP1 constraint between -Z and y'
% a6.bodyJ = 1;
% a6.bodyI = 2;
% a6.aBarJ = [0 0 -1]';
% a6.aBarI = [0 1 0]';
% a6.ft = @(t)cos((pi*cos(2*t))/4 + pi/2);
% a6.ftDot = @(t)((pi*sin(2*t)*sin((pi*cos(2*t))/4 + pi/2))/2);
% a6.ftDDot = @(t)(pi*cos(2*t)*sin((pi*cos(2*t))/4 + pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 + pi/2))/4);
% a6.constraintName = 'DP1 driving constraint';
% isKinematic = 0;
% sys.addBasicConstraint(isKinematic,'dp1',a6);

%% Perform dynamics analysis
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    dynamicsAnalysisTime = toc;
    save('testRevoluteJoint.mat','sys');
else
    load('testRevoluteJoint.mat')
end

disp(['Dynamics Analysis for A8P1 took ' num2str(dynamicsAnalysisTime) ' seconds.'])

%% Display torque at the revolute joint 
% Extract torque for body 2 due to all constraints and time
torque = sys.myBodies{2}.myConstraintTorquesOmegaTotal;
time = sys.myBodies{2}.myTimeTotal;

% Compute theta over time.
theta = pi/4*cos(2*time);

% Extract torque due to DP1 driving constraint.
DP1const = 6;
torqueDriving = torque((3*DP1const-2):3*DP1const,:);

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
plot(time,torqueDriving(3,:));
xlabel('Time (sec)')
ylabel('Torque (N*m)')
axis([0 10 -250 250]);
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