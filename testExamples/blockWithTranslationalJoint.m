%% blockWithTranslationalJoint.m
% File for dynamic simulation of a block with a translational joint. This
% is a simple test case to check implementation of the translational joint.
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 1;
length1 = 0;
isGround1 = 1;
JMatrix1 = eye(3,3);
gravityDirection = '-y';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

% Add body2. Block.
length2 = 0.0173; % The length of the slider does not matter so this is an arbitrary length.
mass2 = 2.0;

% Define inertial properties of the slider
Jxx = 0.0001;
Jyy = 0.0001;
Jzz = 0.0001;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'block', isGround2, mass2, length2, JMatrix2, gravityDirection);

%% Define translational joint
%   body1 = first body in joint
%   body2 = second body in joint
%   pointOnBody1 = 3D location of a point along the translational axis of body1
%   pointOnBody2 = 3D location of a point along the translational axis of body2
%   vector1OnBody1 = 1st vector to define the plane of the
%   cylindrical joint on body1
%   vector2OnBody1 = 2nd vector to define the plane of the
%   cylindrical joint on body1. This is orthogonal to
%   vector1OnBody1.
%   vector1OnBody2 = Vector on body2 that is orthogonal to the
%   plane on body1.
%   vector2OnBody2 = Vector on body2 that is parallel to the
%   plane on body1 and orthogonal to vector1OnBody1.
a.body1 = 1;
a.body2 = 2;
a.pointOnBody1 = [0 0 0]';
a.pointOnBody2 = [1 0 0]';
theta = pi/4;
a.vector1OnBody1 = [sin(theta) cos(theta) 0]';
a.vector2OnBody1 = [0 0 1]';
a.vector1OnBody2 = [1 0 0]';
a.vector2OnBody2 = [0 0 1]';
a.constraintName = 'Translational Joint';

sys.addJoint('translational',a);

%% Set initial conditions of each body
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
% Initial position
r1Initial = zeros(3,1);
r2Initial = zeros(3,1);
rInitial = [r1Initial; r2Initial];

% Initial orientation
p1 = [1 0 0 0]'; 

A2 = [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    0 0 1];
p2 = simEngine3DUtilities.A2p(A2);
pInitial = [p1; p2];
assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( [], rDotInitial, [], [], pDotInitial, [], t);

%% Plot starting configuration of system
sys.plot(1);
% view([90 0])
% saveas(gcf,'A8P1_MechanismInitialPosition.png');

%% Perform dynamics analysis
if 1
    timeStart = 0;
    timeEnd = 5;
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

%% Display position, velocity, and acceleration of the block
% Extract torque for body 2 due to all constraints and time
blockPosition = sys.myBodies{2}.myRTotal;
blockVelocity = sys.myBodies{2}.myRDotTotal;
blockAccel = sys.myBodies{2}.myRDDotTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,blockPosition(1,:))
plot(time,blockPosition(2,:))
plot(time,blockPosition(3,:));
xlabel('Time (sec)')
ylabel('Position (m)')
% h2.Color = 'g';
legend('x','y','z')
title('Block Position')

figure
hold on
plot(time,blockVelocity(1,:))
plot(time,blockVelocity(2,:))
plot(time,blockVelocity(3,:));
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
% h2.Color = 'g';
legend('x','y','z')
title('Block Velocity')

figure
hold on
plot(time,blockAccel(1,:))
plot(time,blockAccel(2,:))
plot(time,blockAccel(3,:));
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
% h2.Color = 'g';
legend('x','y','z')
title('Block Acceleration')

%% Animate the system
plot.animateSystem(sys,[89 -90])
