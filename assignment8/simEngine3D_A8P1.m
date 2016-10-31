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
rInitial = [0 0;
    0 sqrt(2);
    0 -sqrt(2)];

% Initial orientation
p1 = [1 0 0 0]'; 
s2 = sqrt(2)/2;
A2 = [0 0 1;
    s2 s2 0;
    -s2 s2 0];
p2 = simEngine3DUtilities.A2p(A2);
pInitial = [p1 p2];

% Initial velocities. Assume system starts from rest.
rDotInitial = zeros(3,2);
pDotInitial = zeros(4,2);

t = 0;
sys.updateSystemState( rInitial, rDotInitial, [], pInitial, pDotInitial, [], t);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])

%% Build revolute joint from basic constraints
% In future development, include a joint class that automatically creates
% joints from the basic constraints. I will use the definition of a
% revolute joint in an ADAMS input file to guide the user inputs for my
% revolute joint.

% CD constraint for x-value of point Q
a1.bodyI = 1;
a1.bodyJ = 2;
a1.coordVec = [1 0 0]';
a1.sBarJQ = [-2 0 0]';
a1.sBarIP = [0 0 0]';
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
a4.bodyI = 1;
a4.bodyJ = 2;
a4.aBarI = [0 1 0]';
a4.aBarJ = [0 0 1]';
a4.ft = @(t)0;
a4.ftDot = @(t)0;
a4.ftDDot = @(t)0;
a4.constraintName = 'DP1 b/w Y and zPrime';
sys.addBasicConstraint(isKinematic,'dp1',a4);

% DP1 constraint between Z and z'
a5.bodyI = 1;
a5.bodyJ = 2;
a5.aBarI = [0 0 1]';
a5.aBarJ = [0 0 1]';
a5.ft = @(t)0;
a5.ftDot = @(t)0;
a5.ftDDot = @(t)0;
a5.constraintName = 'DP1 b/w Z and zPrime';
sys.addBasicConstraint(isKinematic,'dp1',a5);

%% Perform dynamics analysis
if 1
    timeStart = 0;
    timeEnd = 1;
    timeStep = 10^-3; 
    displayFlag = 0;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, displayFlag);
    save('multibodySystem_A8P1.mat','sys');
else
    load('multibodySystem_A8P1.mat')
end




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
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ')
savefig('A7P1_TorqueVsTime.png')

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
[ax, h1, h2] = plotyy(time,torqueDriving(3,:),time,theta);
xlabel('Time (sec)')
ylabel(ax(1),'Torque (N*m)')
ylabel(ax(2),'Theta (rad)')
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ','Theta')
savefig('A7P1_TorqueAndThetaVsTime.png')
