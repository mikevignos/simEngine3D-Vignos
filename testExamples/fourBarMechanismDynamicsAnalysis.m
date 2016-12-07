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
if exist('fourBarMechanismDynamicsAnalysis.mat')
    load('fourBarMechanismDynamicsAnalysis.mat')
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
    knownInitialOmegaBar = [1; 0; 0];
    knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);
    
    sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);
%     sys.computeAndSetInitialVelocities([], [], []);
%     save('fourBarMechanismDynamicsAnalysis.mat','sys');
end

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 2.5;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    method = 'quasiNewton';
    tic;
%         sys.inverseDynamicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
%         sys.kinematicsAnalysis(timeStart, timeEnd, timeStep, displayFlag);
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag);
    analysisTime = toc;
    save('uniqueFourBarMechanism.mat','sys');
else
    load('uniqueFourBarMechanism.mat')
end

disp(['Inverse Dynamics Analysis for uniqueFourBarMechanism took ' num2str(analysisTime) ' seconds.'])

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

%% Plot the position, velocity, and acceleration of the rocker vs time
rockerPosition = sys.myBodies{4}.myRTotal;
rockerVelocity = sys.myBodies{4}.myRDotTotal;
rockerAccel = sys.myBodies{4}.myRDDotTotal;


figure
hold on
plot(time,rockerPosition(1,:))
plot(time,rockerPosition(2,:))
plot(time,rockerPosition(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Position (m)')
title('Rocker Position')
hold off


figure
hold on
plot(time,rockerVelocity(1,:))
plot(time,rockerVelocity(2,:))
plot(time,rockerVelocity(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
title('Rocker Velocity')
hold off


figure
hold on
plot(time,rockerAccel(1,:))
plot(time,rockerAccel(2,:))
plot(time,rockerAccel(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Rocker Acceleration')
hold off

%% Display torque at the revolute joint 
% Extract torque for body 2 due to all constraints and time
torque = sys.myBodies{2}.myConstraintTorquesOmegaTotal;
time = sys.myBodies{2}.myTimeTotal;

% Compute theta over time.
theta = 2*pi*time + pi/2;

% Extract torque due to DP1 driving constraint.
DP1const = 18;
torqueDriving = torque((3*DP1const-2):3*DP1const,:);

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
plot(time,torqueDriving(3,:));
xlabel('Time (sec)')
ylabel('Torque (N*m)')
% axis([0 1 -1 1]);
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ')
title('Torque Due to DP1 Driving Constraint in Pendulum Reference Frame')
saveas(gcf,'A8P1_TorqueVsTime.png')

figure
hold on
plot(time,torqueDriving(1,:))
plot(time,torqueDriving(2,:))
[ax, h1, h2] = plotyy(time,torqueDriving(3,:),time,theta);
ax(1).YLim = [-250 250];
xlabel('Time (sec)')
ylabel(ax(1),'Torque (N*m)')
ylabel(ax(2),'Theta (rad)')
% h2.Color = 'g';
legend('TorqueX','TorqueY','TorqueZ','Theta')
saveas(gcf,'A8P1_TorqueAndThetaVsTime.png')

%% Display the local reaction x force at the origin of the wheel
force = sys.myBodies{2}.myConstraintForcesTotal;
time = sys.myBodies{2}.myTimeTotal;

% Compute the total force due to the revolute joint constraints
revJointConsts = 1:5;
reactionForceTotal = zeros(3,length(time));
for iC = 1:length(revJointConsts);
    reactionForceTotal = reactionForceTotal + force((3*iC-2):3*iC,:);
end

figure
hold on
plot(time, -1*reactionForceTotal(1,:))
xlabel('Time (sec)')
ylabel('Reaction Force (N)')
title('Reaction Force in X-direction due to revolute joint b/w wheel and ground')
axis([0 2.5 -200 40])
hold off
    
