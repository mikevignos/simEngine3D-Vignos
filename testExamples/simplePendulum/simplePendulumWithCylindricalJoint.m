%% simplePendulumWithCylindricalJoint.m
% Taken from: http://lim.ii.udc.es/mbsbenchmark/dist/A01/A01_specification.xml.
% Used to validate implementation of a cylindrical joint.
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
length2 = 0.0; % meters
mass2 = 1.0; % kg

% Define inertial properties of the mass
Jxx = 0.0;
Jyy = 0.0;
Jzz = 0.0;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'point', isGround2, mass2, length2, JMatrix2, gravityDirection);

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

%% Define cylindrical joint between ground and point mass
%   body1 = first body in joint
%   body2 = second body in joint
%   pointOnBody1 = 3D location of a point along the translational axis of body1
%   pointOnBody2 = 3D location of a point along the translational axis of body2
%   vector1OnBody1 = 1st vector to define the plane of the
%   cylindrical joint on body1
%   vector2OnBody1 = 2nd vector to define the plane of the
%   cylindrical joint on body1.
%   vectorOnBody2 = Vector on body2 that is orthogonal to the
%   plane on body1.
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 1]';
a1.pointOnBody2 = [1 0 0]';
a1.vector1OnBody1 = [1 0 0]';
a1.vector2OnBody1 = [0 1 0]';
a1.vectorOnBody2 = [0 0 1]';
a1.constraintName = 'Cylindrical Joint b/w Ground and Mass';

sys.addJoint('cylindrical',a1);

%% Add force to point mass in z-direction
sys.addForce(2,[0 0 1]', [0 0 0]', 'Force');

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 20.0;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, 0);
    analysisTime = toc;
    save('simplePendulumWithCylindricalJoint.mat','sys');
else
    load('simplePendulumWithCylindricalJoint.mat')
end

disp(['Dynamics Analysis for simplePendulumWithCylindricalJoint took ' num2str(analysisTime) ' seconds.'])

%% Animate system
plot.animateSystem(sys, [0 90]);

%% Plot the position of point mass versus time
massPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,massPosition(1,:));
plot(time,massPosition(2,:));
legend('x','y')
hold off


figure
hold on
plot(time,massPosition(3,:));
title('Z-position')
hold off

%% Plot mass acceleration
massPosition = sys.myBodies{2}.myRTotal;
massVel = sys.myBodies{2}.myRDotTotal;
massAccel = sys.myBodies{2}.myRDDotTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
subplot(3,1,1)
plot(time,massPosition(3,:),'LineWidth',2);
ylabel('Position (m)','FontSize',14)
set(gca,'FontSize',10)
title('Position, Velocity, and Acceleration of Mass in Z-Direction','FontSize',16)

subplot(3,1,2)
plot(time,massVel(3,:),'LineWidth',2);
ylabel('Velocity (m/s)','FontSize',14)
set(gca,'FontSize',10)

subplot(3,1,3)
plot(time,massAccel(3,:),'LineWidth',2);
ylabel('Acceleration (m/s^2)','FontSize',14)
xlabel('Time (sec)','FontSize',14)
set(gca,'FontSize',10)
hold off


%% Compare the results to the benchmark results
benchmark = dlmread('A01_solution_data.txt','\t',1,0);
benchmarkTime = benchmark(:,1);
position = benchmark(:,2:3);
massPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,massPosition(1,:),'r');
plot(time,massPosition(2,:),'b');
plot(benchmarkTime(1:41), position(1:41,1),'rs')
plot(benchmarkTime(1:41), position(1:41,2),'bs')
xlabel('Time (sec','FontSize',16)
ylabel('Position (m)','FontSize',16);
legend('Simulation X','Simulation Y','Benchmark X','Benchmark Y')
set(gca,'FontSize',12)
