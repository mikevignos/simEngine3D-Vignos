%% fourBarMechanismDynamicsAnalysis.m
% Taken from http://lim.ii.udc.es/mbsbenchmark/dist/A02/A02_specification.xml
% Used to test the combinaton of a revolute, universal, and spherical joint.  .
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

% Add body2. The is the first link.
length2 = 1.0; % meters
mass2 = 1.0; % kg

% Define inertial properties of the first link
Jxx = mass2*(length2)^2/12;
Jyy = 0.0;
Jzz = mass2*(length2)^2/12;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2, gravityDirection);

% Add body3. The is the second link (horizontal link).
length3 = 1.0; % meters
mass3 = 1.0; % kg

% Define inertial properties of the first link
Jxx = mass3*(length3)^2/12;
Jyy = 0.0;
Jzz = mass3*(length3)^2/12;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3, gravityDirection);

% Add body4. The is the third link.
length4 = 1.0; % meters
mass4 = 1.0; % kg

% Define inertial properties of the first link
Jxx = mass4*(length4)^2/12;
Jyy = 0.0;
Jzz = mass4*(length4)^2/12;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'bar', isGround4, mass4, length4, JMatrix4, gravityDirection);

%% Define revolute joint between first link and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 0]';
a1.pointOnBody2 = [0 -0.5 0]';
a1.vector1OnBody1 = [1 0 0]';
a1.vector2OnBody1 = [0 1 0]';
a1.vectorOnBody2 = [0 0 1]';
a1.constraintName = 'Revolute joint b/w ground and first link';

sys.addJoint('revolute',a1);

%% Define revolute joint between first link and second link
a2.body1 = 2;
a2.body2 = 3;
a2.pointOnBody1 = [0 0.5 0]';
a2.pointOnBody2 = [0 0.5 0]';
a2.vector1OnBody1 = [1 0 0]';
a2.vector2OnBody1 = [0 1 0]';
a2.vectorOnBody2 = [0 0 1]';
a2.constraintName = 'Revolute joint b/w first and second links';

sys.addJoint('revolute',a2);

%% Define revolute joint between second link and third link
a3.body1 = 3;
a3.body2 = 4;
a3.pointOnBody1 = [0 -0.5 0]';
a3.pointOnBody2 = [0 0.5 0]';
a3.vector1OnBody1 = [1 0 0]';
a3.vector2OnBody1 = [0 1 0]';
a3.vectorOnBody2 = [0 0 1]';
a3.constraintName = 'Revolute joint b/w second and third links';

sys.addJoint('revolute',a3);

%% Define revolute joint between third link and ground
a4.body1 = 1;
a4.body2 = 4;
a4.pointOnBody1 = [1 0 0]';
a4.pointOnBody2 = [0 -0.5 0]';
a4.vector1OnBody1 = [1 0 0]';
a4.vector2OnBody1 = [0 1 0]';
a4.vectorOnBody2 = [0 0 1]';
a4.constraintName = 'Revolute joint b/w third link and ground';

sys.addJoint('revolute',a4);

%% Set initial conditions of each body
if exist('fourBarMechanismDynamicsAnalysisOverConstrained.mat')
    load('fourBarMechanismDynamicsAnalysisOverConstrained.mat')
else
    r1Initial = zeros(3,1); % Ground
    r2Initial = [0.0; 0.5; 0.0]; % Link1
    r3Initial = [0.5; 1.0; 0.0]; % Link2
    r4Initial = [1.0; 0.5; 0.0]; % Link3
    rInitial = [r1Initial; r2Initial; r3Initial; r4Initial];
    
    % Initial orientation
    % Ground
    p1 = [1.0 0.0 0.0 0.0]';
    
    % Link1
    p2 = [1.0 0.0 0.0 0.0]';
    
    % Link2
    A3 = [0 -1 0;
        1 0 0;
        0 0 1];
    p3 = simEngine3DUtilities.A2p(A3);
    
    % Link3
    p4 =  [1 0 0 0]';
    
    pInitial = [p1; p2; p3; p4];
    
    assemblyAnalysisFlag = 1;
    sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);
    
    % Plot starting configuration of system
    sys.plot(1);
    
    % Initial velocities
    known = 2;
    knownInitialRDot = [0.5 0 0]';
    knownInitialOmegaBar = [0; 0; -1];
    knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);
    
    sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);
%     sys.computeAndSetInitialVelocities([], [], []);
%     save('fourBarMechanismDynamicsAnalysis.mat','sys');
end

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 10.0;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
%         sys.inverseDynamicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
%         sys.kinematicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    analysisTime = toc;
    save('fourBarMechanismDynamicsAnalysisOverConstrained.mat','sys');
else
    load('fourBarMechanismDynamicsAnalysisOverConstrained.mat')
end

disp(['Inverse Dynamics Analysis for fourBarMechanismDynamicsAnalysisOverConstrained took ' num2str(analysisTime) ' seconds.'])

%% Plot the position of point B0 versus time
linkOrientation = sys.myBodies{2}.myPTotal;
linkPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;
pointBPosition = zeros(3,length(time));


sB = [0 0.5 0]';
for iT = 1:length(time)
    A = simEngine3DUtilities.p2A(linkOrientation(:,iT));
    pointBPosition(:,iT) = linkPosition(:,iT) + A*sB;
end

figure
hold on
plot(time,pointBPosition(1,:),'r');
plot(time,pointBPosition(2,:),'b');
plot(time,pointBPosition(3,:),'g');
legend('x','y','z');
xlabel('Time (sec)');
ylabel('Displacement (m)');
hold off
