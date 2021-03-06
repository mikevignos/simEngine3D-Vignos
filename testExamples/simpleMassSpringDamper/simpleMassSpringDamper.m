%% pendulumWithRevoluteJoint.m
% File for dynamic simulation of simple pendulum.
% Compared to solution of second order ODE solved using ODE45 in MATLAB
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 0;
length1 = 0;
isGround1 = 1;
JMatrix1 = zeros(3,3);
gravityDirection = '-y';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

% Add body2. This is just a point mass.
length2 = 0; % meters
mass2 = 1; %kg.

% Compute JMatrix.
Jxx = 0.0001;
Jyy = 0.0001;
Jzz = 0.0001;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'point', isGround2, mass2, length2, JMatrix2, gravityDirection);

%% Add TSDA between Ground and Mass
a.bodyI = 1;
a.bodyJ = 2;
a.sBarIP = [0 0 0]';
a.sBarJQ = [0 0 0]';
a.stiffness = 100;
a.restingLength = 0;
a.dampingCoefficient = 2;
a.actuatorFunction = @(lij, lijDot, t)0;
a.name = 'Spring-Damper Element';
sys.addTSDA(a);


%% Define important points and vectors in the system
% sys.addPoint(1, [0 0 0]', 'O');
% sys.addPoint(2, [-2 0 0]', 'Q');

%% Set initial conditions of each body
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
% Initial position
r1Initial = zeros(3,1);
r2Initial = [0.02; 0; 0];
rInitial = [r1Initial; r2Initial];

% Initial orientation
p1 = [1 0 0 0]'; 
p2 = [1 0 0 0]';
pInitial = [p1; p2];

assemblyAnalysisFlag = 1;
sys.setInitialPose(rInitial, pInitial, assemblyAnalysisFlag);

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( [], rDotInitial, [], [], pDotInitial, [], t);

%% Plot starting configuration of system
sys.plot(1);
% view([90 0])

%% Define translation joint
a.body1 = 1;
a.body2 = 2;
a.pointOnBody1 = [0 0 0]';
a.pointOnBody2 = [3 0 0]';
a.vector1OnBody1 = [0 1 0]';
a.vector2OnBody1 = [0 0 1]';
a.vector1OnBody2 = [1 0 0]';
a.vector2OnBody2 = [0 0 1]';
a.constraintName = 'Translational Joint';

sys.addJoint('translational',a);


%% Perform dynamics analysis
if 1
    timeStart = 0;
    timeEnd = 4;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    velocityConstFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velocityConstFlag);
    dynamicsAnalysisTime = toc;
    save('testRevoluteJoint.mat','sys');
else
    load('testRevoluteJoint.mat')
end

disp(['Dynamics Analysis for A8P1 took ' num2str(dynamicsAnalysisTime) ' seconds.'])

%% Animate system
plot.animateSystem(sys, [0 90]);

%% Plot position of mass versus time
position = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,position(1,:));
plot(time,position(2,:));
plot(time,position(3,:));
xlabel('Time (sec)')
ylabel('Position (m)')
legend('x','y','z')
title('Position of mass')
hold off

%% Compare to analytical solution
% Analytical solution
tspan=[0 4];
y0=[0.02;0];
[t,y]=ode45('unforced1',tspan,y0);


position = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;

figure
hold on
plot(time,position(1,:),'LineWidth',2);
plot(t(1:2:end),y((1:2:end),1),'ks');
xlabel('Time (sec)','FontSize',16)
ylabel('Position (m)','FontSize',16)
legend('Simulation','Analytical Solution')
% title('Position of mass','FontSize',16)
set(gca,'FontSize',12);
hold off



