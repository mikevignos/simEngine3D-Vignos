%% simplePendulum.m
% Taken from: http://lim.ii.udc.es/mbsbenchmark/dist/A01/A01_specification.xml.
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

% Add body2. The is the mass.
length2 = 0.0; % meters. This is actually the radius of the wheel.
mass2 = 1.0; % kg

% Define inertial properties of the mass
Jxx = 0.1;
Jyy = 0.1;
Jzz = 0.1;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'point', isGround2, mass2, length2, JMatrix2, gravityDirection);

%% Define revolute joint between ground and point mass
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 0]';
a1.pointOnBody2 = [1 0 0]';
a1.vector1OnBody1 = [1 0 0]';
a1.vector2OnBody1 = [0 1 0]';
a1.vectorOnBody2 = [0 0 1]';
a1.constraintName = 'Revolute Joint b/w Ground and Mass';

sys.addJoint('revolute',a1);

%% Set initial conditions of each body
r1Initial = zeros(3,1); % Ground
r2Initial = [-1; 0; 0]; % Point Mass

rInitial = [r1Initial; r2Initial];

% Initial orientation
% Ground
p1Initial = [1.0 0.0 0.0 0.0]'; 

% Point Mass
p2Initial = [1.0 0.0 0.0 0.0]';

pInitial = [p1Initial; p2Initial];

% Initial velocities. Starting from rest
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Plot starting configuration of system
sys.plot(1);

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 10.0;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    velocityConstraintViolationFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velocityConstraintViolationFlag);
    analysisTime = toc;
    save('simplePendulum.mat','sys');
else
    load('simplePendulum.mat')
end

disp(['Dynamics Analysis for simplePendulum took ' num2str(analysisTime) ' seconds.'])

%% Animate system
plot.animateSystem(sys, [0 90]);

%% Plot the position of point mass versus time
massPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,massPosition(1,:));
plot(time,massPosition(2,:));
plot(time,massPosition(3,:));
legend('x','y','z')
hold off