%% simEngine3D_A6P1.m
% File for testing implementation of DP2 and D constraints
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1
mass1 = 10;
length1 = 2;
isGround1 = 0;
sys.addBody(1, 'bar', isGround1, mass1,length1);

% Add body2
mass2 = 10;
length2 = 2;
isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2,length2);

%% Define important points and vectors in the system
sys.addPoint(1, [1 0 0]', 'P');
sys.addPoint(2, [1 0 0]', 'Q');
sys.addVector(1, [1 0 0]', 'aBar1');

%% Update initial conditions of each body
% In a simulation, the state of each body will be updated at each time
% step. Each column represents a different body.
rInitial = [1 1;
    0.5 0;
    0 0];
rDotInitial = zeros(3,2);
pInitial = [1 1
    0 0
    0 0
    0 0];
pDotInitial = zeros(4,2);
t = 0;
sys.updateSystemState( rInitial, rDotInitial, pInitial, pDotInitial, t);

%% Define constraints between bodies
% D constraint
a.bodyI = 1;
a.bodyJ = 2;
a.sBarIP = [1 0 0]';
a.sBarJQ = [1 0 0]';
a.ft = @(t) 0.5^2;
a.ftDot  = @(t) 0;
a.ftDDot  = @(t) 0;
a.constraintName = 'Bar1 and Bar2 D Constraint';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'d',a);

% DP2 constraint
a.bodyI = 1;
a.bodyJ = 2;
a.aBarI = [1 0 0]';
a.sBarIP = [1 0 0]';
a.sBarJQ = [1 0 0]';
a.ft = @(t) 0;
a.ftDot  = @(t) 0;
a.ftDDot  = @(t) 0;
a.constraintName = 'P and Q distance constraint';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'dp2',a);

%% Plot the current position of the system
sys.plot(1);

%% Extract values of constraints
% For testing purposes
% Compute all quantities for the D constraint
sys.computeConstraintProperties(1, t, 1, 1, 1, 1, 1);

% Compute all quantities for the DP2 constraint
sys.computeConstraintProperties(2, t, 1, 1, 1, 1, 1);

%% Future development
% Define joints between bodies
% The idea is to have the joints class call the constraints classes and build
% each joint from basic constraints