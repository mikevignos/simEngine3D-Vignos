%% sliderCrankMechanism.m
% Taken from ch. 12.2 section Edward J. Haug: Computer Aided Kinematics and Dynamics of
% Mechanical Systems (Allyn and Bacon, 1989)
% Used to test the combinaton of a spherical, revolute-cylindrical,
% revolute, and translation joint.
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define bodies and properties of bodies in system
% Add body1. This is the ground.
mass1 = 1;
length1 = 0;
isGround1 = 1;
JMatrix1 = eye(3,3);
gravityDirection = '-z';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

% Add body2. The is the crank
length2 = 0.08; % meters.
mass2 = 0.12; % kg

% Define inertial properties of the crank
Jxx = 0.0001;
Jyy = 0.00001;
Jzz = 0.0001;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2, gravityDirection);

% Add body3 to the system. This is the connector between the slider and the
% crank.
length3 = 0.3; % meters
mass3 = 0.5; % kg

% Define inertial properties of the connecting rod
Jxx = 0.004;
Jyy = 0.0004;
Jzz = 0.004;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3, gravityDirection);

% Add the slider (body4) to the mechanism.
length4 = 0.0173;
mass4 = 2.0;

% Define inertial properties of the slider
Jxx = 0.0001;
Jyy = 0.0001;
Jzz = 0.0001;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'block', isGround4, mass4, length4, JMatrix4, gravityDirection);

%% Define important points and vectors in the system
sys.addPoint(2, [0 0.08 0]', 'B');
sys.addPoint(2, [0 0 0]', 'revJoint');
sys.addPoint(3, [0.15 0 0]', 'C');
sys.addPoint(3, [-0.15 0 0]', 'B');
sys.addPoint(4, [0.0173/2 0 0]','D');
sys.addPoint(4, [-0.0173/2 0 0]','-D');

%% Define revolute joint between crank and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0.1 0.12]';
a1.pointOnBody2 = [0 0 0]';
a1.vector1OnBody1 = [0 1 0]';
a1.vector2OnBody1 = [0 0 1]';
a1.vectorOnBody2 = [1 0 0]';
a1.constraintName = 'Revolute Joint b/w Ground and Crank';

sys.addJoint('revolute',a1);

%% Define spherical joint between connector and crank
a2.body1 = 2;
a2.body2 = 3;
a2.pointOnBody1 = [0.0 0.08 0.0]';
a2.pointOnBody2 = [-0.15 0.0 0.0]';
a2.constraintName = 'Spherical joint b/w connector and crank';

sys.addJoint('spherical',a2);

%% Define revolute-cylindrical composite joint between connector and slider
a3.body1 = 3;
a3.body2 = 4;
a3.pointOnBody1 = [0.15 0 0]';
a3.pointOnBody2 = [1 0 0]';
a3.vectorOnBody1 = [0 1 0]';
a3.vector1OnBody2 = [1 0 0]';
a3.vector2OnBody2 = [0 1 0]';
a3.vector3OnBody2 = [0 0 1]';
a3.constraintName = 'Rev-cylindrical joint b/w connector and slider';

sys.addJoint('revolute-cylindrical',a3);

%% Define translation joint between slider and ground
a4.body1 = 1;
a4.body2 = 4;
a4.pointOnBody1 = [2 0 0]';
a4.pointOnBody2 = [0.1 0 0]';
a4.vector1OnBody1 = [0 1 0]';
a4.vector2OnBody1 = [0 0 1]';
a4.vector1OnBody2 = [1 0 0]';
a4.vector2OnBody2 = [0 0 1]';
a4.constraintName = 'Translation joint b/w slider and ground';

sys.addJoint('translational',a4);

%% Define distance constraint between connector and slider
a5.bodyI = 3;
a5.bodyJ = 4;
a5.sBarIP = [0.15 0 0]';
a5.sBarJQ = [1 0 0]';
a5.ft = @(t)1;
a5.ftDot =  @(t)0;
a5.ftDDot =  @(t)0;
a5.constraintName = 'Distance constraint b/w connector and slider';
isKinematic = 1;
sys.addBasicConstraint(isKinematic,'d',a5);

%% Set initial conditions of each body
r1Initial = zeros(3,1); % Ground
r2Initial = [0.0; 0.1; 0.12]; % Crank
r3Initial = [0.1; 0.05; 0.1]; % Connector
r4Initial = [0.2; 0.0; 0.0]; % Slider.
rInitial = [r1Initial; r2Initial; r3Initial; r4Initial];

% Initial orientation
% Ground
p1 = [1.0 0.0 0.0 0.0]';

% Crank
p2 = [0.7042 0.71 0.0 0.0]';

% Connector (body3)
p3 = [0.8865 -0.21 0.4 -0.1]';

% Slider
p4 =  [1.0 0.0 0.0 0.0]';

pInitial = [p1; p2; p3; p4];

t = 0;
assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);
 
% Initial velocities. Initial velocity known for crank.
known = 2;
knownInitialRDot = [0; 0; 0];
knownInitialOmegaBar = [6; 0; 0];
knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);

sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);

% Use the next command if you are prescribing a driving constraint.
% sys.computeAndSetInitialVelocities([], [], []);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 1;
    timeStep = 10^-3;
    order = 2;
    displayFlag = 1;
    velConstViolFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velConstViolFlag);
    analysisTime = toc;
    save('sliderCrankMechanismDynamicsAnalysisNoDrivingConst.mat','sys');
else
    load('sliderCrankMechanismDynamicsAnalysisNoDrivingConst.mat')
end

disp(['Analysis for sliderCrankMechanismDynamicsAnalysisNoDrivingConst took ' num2str(analysisTime) ' seconds.'])

%% Animate results
plot.animateSystem(sys,[90,0])

%% Plot the position, velocity, and acceleration of the slider vs time
sliderPosition = sys.myBodies{4}.myRTotal;
sliderVelocity = sys.myBodies{4}.myRDotTotal;
sliderAccel = sys.myBodies{4}.myRDDotTotal;
time = sys.myBodies{4}.myTimeTotal;

figure
hold on
plot(time,sliderPosition(1,:),'LineWidth',2)
plot(time,sliderPosition(2,:),'LineWidth',2)
plot(time,sliderPosition(3,:),'LineWidth',2)
legend('x','y','z')
xlabel('Time (sec)','FontSize',16)
ylabel('Position (m)','FontSize',16)
set(gca,'FontSize',12);
axis([0 1 0.15 0.31])
% title('Slider')
hold off


figure
hold on
plot(time,sliderVelocity(1,:),'LineWidth',2)
plot(time,sliderVelocity(2,:),'LineWidth',2)
plot(time,sliderVelocity(3,:),'LineWidth',2)
legend('x','y','z')
xlabel('Time (sec)','FontSize',16)
ylabel('Velocity (m/s)','FontSize',16)
axis([0 1 -1 1])
set(gca,'FontSize',12);
% title('Slider')
hold off


figure
hold on
plot(time,sliderAccel(1,:),'LineWidth',2)
plot(time,sliderAccel(2,:),'LineWidth',2)
plot(time,sliderAccel(3,:),'LineWidth',2)
axis([0 1 -30 40])
legend('x','y','z')
xlabel('Time (sec)','FontSize',16)
ylabel('Acceleration (m/s^2)','FontSize',16)
set(gca,'FontSize',12);
% title('Slider')
hold off

%% Plot the position of point B versus time
crankOrientation = sys.myBodies{2}.myPTotal;
crankPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{2}.myTimeTotal;
pointBPosition = zeros(3,length(time));
sB = [0 0.08 0]';
for iT = 1:length(time)
    A = simEngine3DUtilities.p2A(crankOrientation(:,iT));
    pointBPosition(:,iT) = crankPosition(:,iT) + A*sB;
end

figure
hold on
plot(time,pointBPosition(1,:));
plot(time,pointBPosition(2,:));
plot(time,pointBPosition(3,:));
legend('x','y','z')
title('Position of tip of crank')
hold off