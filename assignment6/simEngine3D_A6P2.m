%% simEngine3D_A6P2.m
% File for testing implementation of DP2 and D constraints
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 0;
length1 = 0;
isGround1 = 1;
sys.addBody(1, 'ground', isGround1, mass1,length1);

% Add body2
mass2 = 10;
length2 = 4;
isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2,length2);

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
t = 0;
sys.updateSystemState( rInitial, [], [], pInitial, [], [], t);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])

%% Build revolute joint from basic constraints
% In future development, include a joint class that automatically creates
% joints from the basic constraints

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
a6.bodyJ = 1;
a6.bodyI = 2;
a6.aBarJ = [0 0 -1]';
a6.aBarI = [1 0 0]';
a6.ft = @(t)cos((pi*cos(2*t))/4);
a6.ftDot = @(t)((pi*sin(2*t)*sin((pi*cos(2*t))/4))/2);
a6.ftDDot = @(t)(pi*cos(2*t)*sin((pi*cos(2*t))/4) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4))/4);
a6.constraintName = 'DP1 driving constraint';
isKinematic = 0;
sys.addBasicConstraint(isKinematic,'dp1',a6);

%% Compute phiFull at t = 0
sys.computePhiFull();
phiFull = sys.myPhiFull;

%% Compute Jacobian, Nu, and qDot
sys.computeQDot(0);
phiFullJacobian = sys.myPhiFullJacobian;
nuTotal = sys.myNu;
rDot = sys.myRDot;
pDot = sys.myPDot;

%% Compute gamma and qDDot
sys.computeQDDot(0);
gammaTotal = sys.myGamma;
rDDot = sys.myRDDot;
pDDot = sys.myPDDot;

%% Display desired parameters
disp('phiFull = ')
disp(phiFull)
disp(' ')
disp('Jacobian of phiFull = ')
disp(phiFullJacobian)
disp(' ')
disp('Nu = ')
disp(nuTotal)
disp(' ')
disp('Gamma = ')
disp(gammaTotal)

