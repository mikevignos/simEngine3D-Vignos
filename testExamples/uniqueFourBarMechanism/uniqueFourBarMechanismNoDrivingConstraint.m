%% sliderCrankMechanism.m
% Taken from ch. 12.3 section Edward J. Haug: Computer Aided Kinematics and Dynamics of
% Mechanical Systems (Allyn and Bacon, 1989)
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
gravityDirection = '-z';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

% Add body2. The is the wheel.
length2 = 2.0; % meters. This is actually the radius of the wheel.
mass2 = 2.0; % kg

% Define inertial properties of the wheel
Jxx = 4.0;
Jyy = 2.0;
Jzz = 0.0;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'wheel', isGround2, mass2, length2, JMatrix2, gravityDirection);

% Add body3 to the system. This is the first link.
length3 = 12.2; % meters
mass3 = 1.0; % kg

% Define inertial properties of the connecting rod
Jxx = 12.4;
Jyy = 0.01;
Jzz = 0.0;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3, gravityDirection);

% Add body4 to the system. This is second link.
length4 = 7.4;
mass4 = 1.0;

% Define inertial properties of body4
Jxx = 4.54;
Jyy = 0.01;
Jzz = 0.0;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'bar', isGround4, mass4, length4, JMatrix4, gravityDirection);

%% Define points on bodies
sys.addPoint(2,[0 0 0]','A');
sys.addPoint(2,[0 0 2]','B');
sys.addPoint(3,[0 6.1 0]','B');
sys.addPoint(3,[0 -6.1 0]','C');
sys.addPoint(4,[0 -3.7 0]','C');
sys.addPoint(4,[0 3.7 0]','D');

%% Define revolute joint between wheel and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 0]';
a1.pointOnBody2 = [0 0 0]';
a1.vector1OnBody1 = [0 1 0]';
a1.vector2OnBody1 = [0 0 1]';
a1.vectorOnBody2 = [1 0 0]';
a1.constraintName = 'Revolute Joint b/w Ground and Wheel';

sys.addJoint('revolute',a1);

%% Define universal joint between wheel and first link
necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} ...
    {'pointOnBody2'} {'vectorOnBody1'} {'vectorOnBody2'}];

a2.body1 = 2;
a2.body2 = 3;
a2.pointOnBody1 = [0.0 0.0 2.0]';
a2.pointOnBody2 = [0.0 6.1 0.0]';
a2.vectorOnBody1 = [0.75 -0.662 0.0]';
a2.vectorOnBody2 = [0 0 1]';
a2.constraintName = 'Universal joint b/w wheel and 1st link';

sys.addJoint('universal',a2);

%% Define spherical joint between the links
a3.body1 = 3;
a3.body2 = 4;
a3.pointOnBody1 = [0.0 -6.1 0.0]';
a3.pointOnBody2 = [0.0 -3.7 0.0]';
a3.constraintName = 'Spherical joint b/w links';

sys.addJoint('spherical',a3);

%% Define revolute joint between second link and ground
a4.body1 = 1;
a4.body2 = 4;
a4.pointOnBody1 = [-4.0 -8.5 0.0]';
a4.pointOnBody2 = [0.0 3.7 0.0]';
a4.vector1OnBody1 = [1 0 0]';
a4.vector2OnBody1 = [0 0 1]';
a4.vectorOnBody2 = [1 0 0]';
a4.constraintName = 'Revolute joint b/w 2nd link and ground';

sys.addJoint('revolute',a4);

%% Add torque
bodyNumber = 2;
torque = [10 0 0]';
torqueName = 'Torque applied to body2';
sys.addTorque(bodyNumber, torque, torqueName);

%% Set initial conditions of each body
r1Initial = zeros(3,1); % Ground
r2Initial = [0.0; 0.0; 0.0]; % Wheel
r3Initial = [-3.75; -4.25; 4.25]; % Link1
r4Initial = [-5.75; -8.5; 3.29]; % Link2
rInitial = [r1Initial; r2Initial; r3Initial; r4Initial];

% Initial orientation
% Ground
p1 = [1.0 0.0 0.0 0.0]';

% Wheel
p2 = [1.0 0.0 0.0 0.0]';

% Link1
p3 = [0.8806 -0.29 -0.27 -0.26]';

% Link2
p4 =  [0.6072 -0.36 0.36 -0.61]';

pInitial = [p1; p2; p3; p4];

assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);

% Initial velocities
% known = 2;
% knownInitialRDot = zeros(3,1);
% knownInitialOmegaBar = [2*pi; 0; 0];
% knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, 2, knownInitialOmegaBar);

% sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);
sys.computeAndSetInitialVelocities([], [], []);

%% Plot starting configuration of system
sys.plot(1);
view([90 0])

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    velConstViolFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velConstViolFlag);
    analysisTime = toc;
    save('uniqueFourBarMechanismNoDrivingConstraint.mat','sys');
else
    load('uniqueFourBarMechanismNoDrivingConstraint.mat')
end

disp(['Inverse Dynamics Analysis for uniqueFourBarMechanismNoDrivingConstraint took ' num2str(analysisTime) ' seconds.'])

%% Animate results
plot.animateSystem(sys,[90,0])

%% Plot the position of point B versus time
wheelOrientation = sys.myBodies{2}.myPTotal;
wheelPosition = sys.myBodies{2}.myRTotal;
time = sys.myBodies{4}.myTimeTotal;
pointBPosition = zeros(3,length(time));


sB = [0 0.0 2.0]';
for iT = 1:length(time)
    A = simEngine3DUtilities.p2A(wheelOrientation(:,iT));
    pointBPosition(:,iT) = wheelPosition(:,iT) + A*sB;
end

figure
hold on
plot(time,pointBPosition(1,:));
plot(time,pointBPosition(2,:));
plot(time,pointBPosition(3,:));
legend('x','y','z')
hold off

%% Plot the Z position and velocity of the rocker vs time
rockerPosition = sys.myBodies{4}.myRTotal;
rockerVelocity = sys.myBodies{4}.myRDotTotal;
rockerAccel = sys.myBodies{4}.myRDDotTotal;

figure
plot(time,rockerPosition(3,:), 'LineWidth',2)
xlabel('Time (sec)','FontSize',16)
ylabel('Z Position (m)','FontSize',16)
title('Rocker Z Position','FontSize',16)
set(gca,'FontSize',12);
axis([timeStart timeEnd 2.5 3.75])


figure
plot(time,rockerVelocity(3,:), 'LineWidth',2)
xlabel('Time (sec)','FontSize',16)
ylabel('Z Velocity (m/s)','FontSize',16)
title('Rocker Z Velocity','FontSize',16)
set(gca,'FontSize',12);
axis([timeStart timeEnd -3 3])
