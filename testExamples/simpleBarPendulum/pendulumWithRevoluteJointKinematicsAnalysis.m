%% pendulumWithRevoluteJoint.m
% File for dynamic simulation of simple pendulum.
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 0;
length1 = 0;
isGround1 = 1;
JMatrix1 = zeros(3,3);
gravityDirection = '-z';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

% Add body2
% Compute mass of the bar.
length2 = 4; % meters
density = 7800; %kg/m^3
area = 0.05 * 0.05; % m^2
volume = length2*area; % m^3
mass2 = density*volume;

% Compute JMatrix.
Jxx = 0;
Jyy = (mass2*length2^2)/12;
Jzz = (mass2*length2^2)/12;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2, gravityDirection);

%% Define important points and vectors in the system
sys.addPoint(1, [0 0 0]', 'O');
sys.addPoint(2, [-2 0 0]', 'Q');

%% Set initial conditions of each body
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
% Initial position
r1Initial = zeros(3,1);
r2Initial = [0; sqrt(2); -sqrt(2)];
rInitial = [r1Initial; r2Initial];

% Initial orientation
p1 = [1 0 0 0]'; 
s2 = sqrt(2)/2;
A2 = [0 0 1;
    s2 s2 0;
    -s2 s2 0];
p2 = simEngine3DUtilities.A2p(A2);
pInitial = [p1; p2];

assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);

% Initial velocities. Assume system starts from rest.
% rDotInitial = zeros(6,1);
% pDotInitial = zeros(8,1);
% 
% t = 0;
% sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

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
a.vector1OnBody1 = [0 1 0]';
a.vector2OnBody1 = [0 0 1]';
a.vectorOnBody2 = [0 0 1]';
a.constraintName = 'Revolute Joint';

sys.addJoint('revolute',a);

%% Add driving constraint to model
% DP1 constraint between -Z and y'
a6.bodyJ = 1;
a6.bodyI = 2;
a6.aBarJ = [0 0 -1]';
a6.aBarI = [0 1 0]';
a6.ft = @(t)cos((pi*cos(2*t))/4 + pi/2);
a6.ftDot = @(t)((pi*sin(2*t)*sin((pi*cos(2*t))/4 + pi/2))/2);
a6.ftDDot = @(t)(pi*cos(2*t)*sin((pi*cos(2*t))/4 + pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 + pi/2))/4);
a6.constraintName = 'DP1 driving constraint';
isKinematic = 0;
sys.addBasicConstraint(isKinematic,'dp1',a6);

%% Perform  analysis
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-3;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
    sys.inverseDynamicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
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