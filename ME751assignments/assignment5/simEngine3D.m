%% simEngine3D.m
% Driver file for dynamics engine
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1
mass1 = 10;
length1 = 1;
isGround1 = 0;
sys.addBody(1, 'bar', isGround1, mass1,length1);

% Add body2
mass2 = 10;
length2 = 1;
isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2,length2);

%% Define important points and vectors in the system
sys.addPoint(1, [1 0 0]', 'P');
sys.addPoint(2, [-0.5 0 0]', 'Q');
sys.addVector(1, [1 0 0]', 'aBar1');
sys.addVector(2, [1 0 0]', 'aBar2');


%% Update initial conditions of each body
% In a simulation, the state of each body will be updated at each time step
rInitial = [-0.5 0.5;
    1 0;
    0 0];
rDotInitial = zeros(3,2);
pInitial = [sqrt(2)/2 1
    0 0
    0 0
    -sqrt(2)/2 0];
pDotInitial = zeros(4,2);
t = 0;
sys.updateSystemState( rInitial, rDotInitial, pInitial, pDotInitial, t);

%% Define constraints between bodies
% DP1 constraint
bodyI = 1;
bodyJ = 2;
aBarI = [1 0 0]';
aBarJ = [1 0 0]';
ft = @(t) 0;
ftDot  = @(t) 0;
ftDDot  = @(t) 0;
constraintName = 'Bar1 and Bar2 Orthogonal Constraint';
sys.addDP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName);

% CD constraint
c = [0 1 0]';
sBarIP = [1 0 0]';
sBarJQ = [-0.5 0 0]';
constraintName = 'P and Q coordinate constraint';
sys.addCDconstraint(bodyI, bodyJ, c, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName);

%% Extract values of constraints
% For testing purposes
% Compute all quantities for the DP1 constraint
sys.computeConstraintProperties(1, t, 1, 1, 1, 1, 1);

% Compute all quantities for the CD constraint
sys.computeConstraintProperties(2, t, 1, 1, 1, 1, 1);


%% Future development
% Define joints between bodies
% The idea is to have the joints class call the constraints classes and build
% each joint from basic constraints