%% simEngine3D_A9P1.m
% File for dynamic simulation of simple pendulum.
clear; close all; clc;
addpath('simEngine3DCode');

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 0;
length1 = 0;
isGround1 = 1;
JMatrix1 = zeros(3,3);
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1);

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
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2);

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

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Plot starting configuration of system
% sys.plot(1);
% view([90 0])
% saveas(gcf,'A9P1_MechanismInitialPosition.png');

%% Build revolute joint from basic constraints
% In future development, include a joint class that automatically creates
% joints from the basic constraints. I will use the definition of a
% revolute joint in an ADAMS input file to guide the user inputs for my
% revolute joint.

% CD constraint for x-value of point Q
a1.bodyJ = 1;
a1.bodyI = 2;
a1.coordVec = [1 0 0]';
a1.sBarIP = [-2 0 0]';
a1.sBarJQ = [0 0 0]';
a1.ft = @(t)0;
a1.ftDot =  @(t)0;
a1.ftDDot =  @(t)0;
a1.constraintName = 'CD constraint on Qx';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a1);

% CD constraint for y-value of point Q
a2 = a1;
a2.coordVec = [0 1 0]';
a2.constraintName = 'CD constraint on Qy';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a2);

% CD constraint for z-value of point Q
a3 = a1;
a3.coordVec = [0 0 1]';
a3.constraintName = 'CD constraint on Qz';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a3);

% DP1 constraint between Y and z'
a4.bodyJ = 1;
a4.bodyI = 2;
a4.aBarJ = [0 1 0]';
a4.aBarI = [0 0 1]';
a4.ft = @(t)0;
a4.ftDot = @(t)0;
a4.ftDDot = @(t)0;
a4.constraintName = 'DP1 b/w Y and zPrime';
sys.addBasicConstraint(isKinematic,'dp1',a4);

% DP1 constraint between Z and z'
a5.bodyJ = 1;
a5.bodyI = 2;
a5.aBarJ = [0 0 1]';
a5.aBarI = [0 0 1]';
a5.ft = @(t)0;
a5.ftDot = @(t)0;
a5.ftDDot = @(t)0;
a5.constraintName = 'DP1 b/w Z and zPrime';
sys.addBasicConstraint(isKinematic,'dp1',a5);

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

%% Perform dynamics analysis with Newton-Raphson method
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 0;
    method = 'newtonRaphson';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    dynamicsAnalysisTime = toc;
    save('multibodySystem_A9P1_NR.mat','sys');
else
    load('multibodySystem_A9P1_NR.mat')
end

disp(['Dynamics Analysis for A9P1 with Newton-Raphson took ' num2str(dynamicsAnalysisTime) ' seconds.'])



%% Extract accelerations of pendulum and number of iterations for Newton-Raphson method
rDDotNR = sys.myBodies{2}.myRDDotTotal;
nIterNR = sys.myIterCountTotal;

%% Reset initial conditions
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

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Perform dynamics analysis with modified Newton method
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 0;
    method = 'modifiedNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    dynamicsAnalysisTime = toc;
    save('multibodySystem_A9P1_modfiedNewton.mat','sys');
else
    load('multibodySystem_A9P1_modfiedNewton.mat')
end

disp(['Dynamics Analysis for A9P1 with modified Newton took ' num2str(dynamicsAnalysisTime) ' seconds.'])

%% Extract accelerations of pendulum and number of iterations for modified Newton method
rDDotModN = sys.myBodies{2}.myRDDotTotal;
nIterModN = sys.myIterCountTotal;

%% Reset initial conditions
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

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(6,1);
pDotInitial = zeros(8,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Perform dynamics analysis with quasi Newton method
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    dynamicsAnalysisTime = toc;
    save('multibodySystem_A9P1_quasiNewton.mat','sys');
else
    load('multibodySystem_A9P1_quasiNewton.mat')
end

disp(['Dynamics Analysis for A9P1 with quasi-Newton took ' num2str(dynamicsAnalysisTime) ' seconds.'])

%% Extract accelerations of pendulum and number of iterations for quasi Newton method
rDDotQuasiN = sys.myBodies{2}.myRDDotTotal;
nIterQuasiN = sys.myIterCountTotal;
time = sys.myTimeTotal;

%% Create plots
rDDotQuasiNvsNR = abs(rDDotQuasiN - rDDotNR);
rDDotModNvsNR = abs(rDDotModN - rDDotNR);

figure
subplot(3,1,1)
plot(time, rDDotQuasiNvsNR(1,:),'r');
ylabel('Difference (m/s^2)')
legend('X');
title('Quasi-Newton vs Newton-Raphson Acceleration');

subplot(3,1,2)
plot(time, rDDotQuasiNvsNR(2,:),'g');
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Y');

subplot(3,1,3)
plot(time, rDDotQuasiNvsNR(3,:),'b');
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Z');
saveas(gcf,'A9P1_AccelerationDiff_QuasiNewtonvsNewtonRaphson.png')

figure
subplot(3,1,1)
plot(time, rDDotModNvsNR(1,:),'r');
title('Modified-Newton vs Newton-Raphson Acceleration')
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('X');

subplot(3,1,2)
plot(time, rDDotModNvsNR(2,:),'g');
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Y');

subplot(3,1,3)
plot(time, rDDotModNvsNR(3,:),'b');
xlabel('Time (sec)')
ylabel('Difference (m/s^2)')
legend('Z');
saveas(gcf,'A9P1_AccelerationDiff_ModNewtonvsNewtonRaphson.png')

figure
hold on
plot(time,nIterNR,'LineWidth',4);
plot(time,nIterModN,'LineWidth',2);
plot(time,nIterQuasiN,'LineWidth',1);
xlabel('Time (sec)')
ylabel('# of Iterations')
title('Number of Iterations to Converge')
legend('Newton-Raphson','Modified-Newton','Quasi-Newton')
axis([0 10 0 5])
hold off
saveas(gcf,'A9P1_IterationCount.png')