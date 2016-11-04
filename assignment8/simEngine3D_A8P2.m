%% simEngine3D_A8P1.m
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
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1);

% Add body2
% Compute mass of the bar.
length2 = 4; % meters
density = 7800; %kg/m^3
area = 0.05 * 0.05; % m^2
volume2 = length2*area; % m^3
mass2 = density*volume2;

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

% Add body3.
% Compute mass of bar
length3 = 2; % meters
density = 7800; %kg/m^3
area = 0.05 * 0.05; % m^2
volume3 = length3*area; % m^3
mass3 = density*volume3;

% Compute polar moment of inertia matrix
Jxx = 0;
Jyy = (mass3*length2^3)/12;
Jzz = (mass3*length2^3)/12;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3);

%% Define important points and vectors in the system
sys.addPoint(1, [0 0 0]', 'O');
sys.addPoint(2, [-2 0 0]', 'Q');
sys.addPoint(2, [2 0 0]', 'P2');
sys.addPoint(3, [-1 0 0]', 'P3');

%% Set initial conditions of each body
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
% Initial position
r1Initial = zeros(3,1);
r2Initial = [0 2 0]';
r3Initial = [0 4 -1]';
rInitial = [r1Initial; r2Initial; r3Initial];

% Initial orientation
p1 = [1 0 0 0]'; 
A2 = [0 0 1;
    1 0 0;
    0 1 0];
p2 = simEngine3DUtilities.A2p(A2);

A3 = [0 0 1;
    0 1 0;
    -1 0 0];
p3 = simEngine3DUtilities.A2p(A3);
pInitial = [p1; p2; p3];

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(9,1);
pDotInitial = zeros(12,1);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])

%% Build revolute joint between ground and body2 from basic constraints
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

% DP1 constraint between Y and z2'
a4.bodyJ = 1;
a4.bodyI = 2;
a4.aBarJ = [0 1 0]';
a4.aBarI = [0 0 1]';
a4.ft = @(t)0;
a4.ftDot = @(t)0;
a4.ftDDot = @(t)0;
a4.constraintName = 'DP1 b/w Y and z2Prime';
sys.addBasicConstraint(isKinematic,'dp1',a4);

% DP1 constraint between Z and z2'
a5.bodyJ = 1;
a5.bodyI = 2;
a5.aBarJ = [0 0 1]';
a5.aBarI = [0 0 1]';
a5.ft = @(t)0;
a5.ftDot = @(t)0;
a5.ftDDot = @(t)0;
a5.constraintName = 'DP1 b/w Z and z2Prime';
sys.addBasicConstraint(isKinematic,'dp1',a5);

%% Build revolute joint between body2 and body3 from basic constraints
% In future development, include a joint class that automatically creates
% joints from the basic constraints. I will use the definition of a
% revolute joint in an ADAMS input file to guide the user inputs for my
% revolute joint.

% CD constraint for x-value of point P
a6.bodyJ = 3;
a6.bodyI = 2;
a6.coordVec = [1 0 0]';
a6.sBarIP = [2 0 0]';
a6.sBarJQ = [-1 0 0]';
a6.ft = @(t)0;
a6.ftDot =  @(t)0;
a6.ftDDot =  @(t)0;
a6.constraintName = 'CD constraint on Px';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a6);

% CD constraint for y-value of point P
a7 = a6;
a7.coordVec = [0 1 0]';
a7.constraintName = 'CD constraint on Py';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a7);

% CD constraint for z-value of point P
a8 = a6;
a8.coordVec = [0 0 1]';
a8.constraintName = 'CD constraint on Pz';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'cd',a8);

% DP1 constraint between x2' and z3'
a9.bodyJ = 3;
a9.bodyI = 2;
a9.aBarJ = [0 0 1]';
a9.aBarI = [1 0 0]';
a9.ft = @(t)0;
a9.ftDot = @(t)0;
a9.ftDDot = @(t)0;
a9.constraintName = 'DP1 b/w x2Prime and z3Prime';
sys.addBasicConstraint(isKinematic,'dp1',a9);

% DP1 constraint between y2' and z3'
a10.bodyJ = 3;
a10.bodyI = 2;
a10.aBarJ = [0 0 1]';
a10.aBarI = [0 1 0]';
a10.ft = @(t)0;
a10.ftDot = @(t)0;
a10.ftDDot = @(t)0;
a10.constraintName = 'DP1 b/w y2Prime and z3Prime';
sys.addBasicConstraint(isKinematic,'dp1',a10);

%% Perform dynamics analysis
if 1
    timeStart = 0;
    timeEnd = 2;
    timeStep = 10^-3;
    order = 1;
    displayFlag = 1;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, displayFlag);
    save('multibodySystem_A8P2.mat','sys');
else
    load('multibodySystem_A8P2.mat')
end

%% Create a plot of the origin of the first pendulum.
time = sys.myTimeTotal;
rOprime1 = sys.myBodies{2}.myRTotal;
rDotOprime1 = sys.myBodies{2}.myRDotTotal;
rDDotOprime1 = sys.myBodies{2}.myRDDotTotal;

figure
subplot(3,1,1)
hold on
plot(time,rOprime1(1,:));
plot(time,rOprime1(2,:));
plot(time,rOprime1(3,:));
xlabel('Time (sec)');
ylabel('Position (m)');
title('Position of origin of 1st pendulum')
legend('X','Y','Z')

subplot(3,1,2)
hold on
plot(time,rDotOprime1(1,:));
plot(time,rDotOprime1(2,:));
plot(time,rDotOprime1(3,:));
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
title('Velocity of origin of 1st pendulum')

subplot(3,1,3)
hold on
plot(time,rDDotOprime1(1,:));
plot(time,rDDotOprime1(2,:));
plot(time,rDDotOprime1(3,:));
xlabel('Time (sec)');
ylabel('Acceleration (m/s^2)');
title('Acceleration of origin of 1st pendulum')

%% Create a plot of the origin of the second pendulum
rOprime2 = sys.myBodies{3}.myRTotal;
rDotOprime2 = sys.myBodies{3}.myRDotTotal;
rDDotOprime2 = sys.myBodies{3}.myRDDotTotal;

figure
subplot(3,1,1)
hold on
plot(time,rOprime2(1,:));
plot(time,rOprime2(2,:));
plot(time,rOprime2(3,:));
xlabel('Time (sec)');
ylabel('Position (m)');
title('Position of origin of 2nd pendulum')
legend('X','Y','Z')

subplot(3,1,2)
hold on
plot(time,rDotOprime2(1,:));
plot(time,rDotOprime2(2,:));
plot(time,rDotOprime2(3,:));
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
title('Velocity of origin of 2nd pendulum')

subplot(3,1,3)
hold on
plot(time,rDDotOprime2(1,:));
plot(time,rDDotOprime2(2,:));
plot(time,rDDotOprime2(3,:));
xlabel('Time (sec)');
ylabel('Acceleration (m/s^2)');
title('Acceleration of origin of 2nd pendulum')

%% Create a plot of the violation of the velocity constraint
velConstViolation = sys.myVelocityConstraintViolationTotal;

% Extract portion that relates to the revolute joint between the pendulums
violation = velConstViolation(6:10,:);

% Compute and plot norm-2
violationNorm = zeros(1,length(time));
for iT = 1:length(time)
    violationNorm(1,iT) = norm(violation(:,iT),2);
end

figure
plot(time, violationNorm)
xlabel('Time')
ylabel('Norm-2 of Vel Const Violation')
title('Norm-2 of Violation of Velocity Constraint for Rev Joint b/w Pendulums')


