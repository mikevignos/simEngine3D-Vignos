%% flyballGovernor.m
% Taken from: http://lim.ii.udc.es/mbsbenchmark/dist/A05/A05_specification.xml
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

% Add body2. The is the rotation shaft.
length2 = 1.0; % meters.
density2 = 3000;
area2 = 0.01*0.01;
volume2 = length2*area2;
mass2 = density2*volume2;

% Define inertial properties of the rotation shaft.
Jxx = mass2*length2^2/12;
Jyy = mass2*length2^2/12;
Jzz = 0;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2, gravityDirection);

% Add body3 to the system. This is the right side rod in the picture.
length3 = 1.0; % meters.
density3 = 3000;
area3 = 0.01*0.01;
volume3 = length3*area3;
mass3 = density3*volume3;

% Define inertial properties of the rotation shaft.
Jxx = mass3*length3^2/12;
Jyy = mass3*length3^2/12;
Jzz = 0;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'bar', isGround3, mass3, length3, JMatrix3, gravityDirection);

% Add body4 to the system. This is the right side rod in the picture.
length4 = 1.0; % meters.
density4 = 3000;
area4 = 0.01*0.01;
volume4 = length4*area4;
mass4 = density4*volume4;

% Define inertial properties of the rotation shaft.
Jxx = mass4*length4^2/12;
Jyy = mass4*length4^2/12;
Jzz = 0;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'bar', isGround4, mass4, length4, JMatrix4, gravityDirection);

% Add body5. This is the base.
length5 = 0.1;
area5 = 0.1*0.1;
volume5 = length5*area5;
density = 3000;
mass5 = volume5*density;

% Define inertial properties of the slider
Jxx = mass5*length5^2/6;
Jyy = mass5*length5^2/6;
Jzz = mass5*length5^2/6;
JMatrix5 = zeros(3,3);
JMatrix5(1,1) = Jxx;
JMatrix5(2,2) = Jyy;
JMatrix5(3,3) = Jzz;

isGround5 = 0;
sys.addBody(5, 'block', isGround5, mass5, length5, JMatrix5, gravityDirection);

%% Add TSDAs to the system
factor = 1e5;
s1.bodyI = 3;
s1.bodyJ = 5;
s1.sBarIP = [0 0 0]';
s1.sBarJQ = [0.05 0 0]';
s1.stiffness = 8e5/factor;
s1.restingLength = 0.5;
s1.dampingCoefficient = 4e4/factor;
s1.actuatorFunction = @(lij, lijDot, t)0;
s1.name = 'Spring-Damper Element b/w body 3 and 5';
sys.addTSDA(s1);
% 
s2.bodyI = 4;
s2.bodyJ = 5;
s2.sBarIP = [0 0 0]';
s2.sBarJQ = [-0.05 0 0]';
s2.stiffness = 8e5/factor;
s2.restingLength = 0.5;
s2.dampingCoefficient = 4e4/factor;
s2.actuatorFunction = @(lij, lijDot, t)0;
s2.name = 'Spring-Damper Element b/w body 4 and 5';
sys.addTSDA(s2);

%% Define important points and vectors in the system
sys.addPoint(2, [0 0 -0.5]', 'base');
sys.addPoint(2, [0 0 0.5]', 'top');
sys.addPoint(3, [0 0 -0.5]', 'base');
sys.addPoint(3, [0 0 0.5]', 'top');
sys.addPoint(4, [0 0 -0.5]', 'base');
sys.addPoint(4, [0 0 0.5]', 'top');
sys.addPoint(5, [0.05 0 0]','connectionForTSDA1');
sys.addPoint(5, [-0.05 0 0]','connectionForTSDA2');

%% Define revolute joint between axis and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 0]';
a1.pointOnBody2 = [0 0 -0.5]';
a1.vector1OnBody1 = [1 0 0]';
a1.vector2OnBody1 = [0 1 0]';
a1.vectorOnBody2 = [0 0 1]';
a1.constraintName = 'Revolute Joint b/w Ground and Axis';

sys.addJoint('revolute',a1);

%% Define translational joint between axis and base
%   body1 = first body in joint
%   body2 = second body in joint
%   pointOnBody1 = 3D location of a point along the translational axis of body1
%   pointOnBody2 = 3D location of a point along the translational axis of body2
%   vector1OnBody1 = 1st vector to define the plane of the
%   translational joint on body1
%   vector2OnBody1 = 2nd vector to define the plane of the
%   translational joint on body1. This is orthogonal to
%   vector1OnBody1.
%   vector1OnBody2 = Vector on body2 that is orthogonal to the
%   plane on body1.
%   vector2OnBody2 = Vector on body2 that is parallel to the
%   plane on body1 and orthogonal to vector1OnBody1.

a2.body1 = 2;
a2.body2 = 5;
a2.pointOnBody1 = [0 0 -0.5]';
a2.pointOnBody2 = [0 0 0]';
a2.vector1OnBody1 = [0 1 0]';
a2.vector2OnBody1 = [1 0 0]';
a2.vector1OnBody2 = [0 0 1]';
a2.vector2OnBody2 = [1 0 0]';
a2.constraintName = 'Translational joint b/w axis and base';

sys.addJoint('translational',a2);

%% Define revolute joint between axis and right rod
%   body1 = first body in joint
%   body2 = second body in joint
%   pointOnBody1 = 3D location of the spherical joint on body1
%   pointOnBody2 = 3D location of the spherical joint on body2
%   vector1OnBody1 = 1st vector to define the plane of the
%   revolute joint on body1
%   vector2OnBody1 = 2nd vector to define the plane of the
%   revolute joiny on body1.
%   vectorOnBody2 = Vector on body2 that is orthogonal to the
%   plane on body1.
a3.body1 = 2;
a3.body2 = 3;
a3.pointOnBody1 = [0.05 0 0.5]';
a3.pointOnBody2 = [0 0 0.5]';
a3.vector1OnBody1 = [1 0 0]';
a3.vector2OnBody1 = [0 0 1]';
a3.vectorOnBody2 = [0 1 0]';
a3.constraintName = 'Revolute Joint b/w right rod and Axis';

sys.addJoint('revolute',a3);
            
%% Define revolute joint between axis and left rod
%   body1 = first body in joint
%   body2 = second body in joint
%   pointOnBody1 = 3D location of the spherical joint on body1
%   pointOnBody2 = 3D location of the spherical joint on body2
%   vector1OnBody1 = 1st vector to define the plane of the
%   revolute joint on body1
%   vector2OnBody1 = 2nd vector to define the plane of the
%   revolute joiny on body1.
%   vectorOnBody2 = Vector on body2 that is orthogonal to the
%   plane on body1.
a4.body1 = 2;
a4.body2 = 4;
a4.pointOnBody1 = [-0.05 0 0.5]';
a4.pointOnBody2 = [0 0 0.5]';
a4.vector1OnBody1 = [1 0 0]';
a4.vector2OnBody1 = [0 0 1]';
a4.vectorOnBody2 = [0 1 0]';
a4.constraintName = 'Revolute Joint b/w left rod and Axis';

sys.addJoint('revolute',a4);

%% Add driving constraint to model
% DP1 constraint between X axis of ground and y axis of body 1.
a5.bodyJ = 1;
a5.bodyI = 2;
a5.aBarJ = [1 0 0]';
a5.aBarI = [0 1 0]';
a5.ft = @(t)cos(2*pi*t + pi/2 + 0.0001);
a5.ftDot = @(t)(-2*pi*sin(pi/2 + 2*pi*t + 0.0001));
a5.ftDDot = @(t)(-4*pi^2*cos(pi/2 + 2*pi*t + 0.0001));
a5.constraintName = 'DP1 driving constraint';
isKinematic = 0;
sys.addBasicConstraint(isKinematic,'dp1',a5);

%% Set initial conditions of each body
% if exist('flyballGovernorSetup.mat')
%     load('flyballGovernorSetup.mat')
% else
alpha = 30;
r1Initial = zeros(3,1); % Ground
r2Initial = [0.0; 0.0; 0.5]; % Axis
r3Initial = [0.05 + 0.5*sind(90-alpha); 0.0; 1 - 0.5*cosd(90-alpha)]; % Right rod
r4Initial = [-0.05 - 0.5*sind(90-alpha); 0.0; 1 - 0.5*cosd(90-alpha)]; % Left rod
r5Initial = [0; 0; 0.5]; % Base
rInitial = [r1Initial; r2Initial; r3Initial; r4Initial; r5Initial];

% Initial orientation
% Ground
p1 = [1.0 0.0 0.0 0.0]';

% Axis
p2 = [1 0 0 0]';

% Right rod
A3 = [cosd(90-alpha) 0 -cosd(alpha);
    0 1 0;
    sind(90-alpha) 0 sind(alpha)];
p3 = simEngine3DUtilities.A2p(A3);

% Left rod
A4 = [cosd(90-alpha) 0 cosd(alpha);
    0 1 0;
    -sind(90-alpha) 0 sind(alpha)];
p4 = simEngine3DUtilities.A2p(A4);

% Base
p5 = [1 0 0 0]';

pInitial = [p1; p2; p3; p4; p5];

t = 0;
assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);
 
% Initial velocities. Initial velocity known for crank.
% known = 2;
% knownInitialRDot = [0; 0; 0];
% knownInitialOmegaBar = [0; 0; 2*pi];
% knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);
% 
% sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);

% Use the next command if you are prescribing a driving constraint.
sys.computeAndSetInitialVelocities([], [], []);


%     save('flyballGovernorSetup.mat','sys');
% end

%% Plot starting configuration of system
sys.plot(1);
view([0 0])

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 1;
    timeStep = 10^-2;
    order = 2;
    displayFlag = 1;
    velocityConstraintFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velocityConstraintFlag);
    analysisTime = toc;
    save('flyballGovernor.mat','sys');
else
    load('flyballGovernor.mat')
end

disp(['Analysis for flyballGovernor took ' num2str(analysisTime) ' seconds.'])

%% Animate results
plot.animateSystem(sys,[0,0])

%% Plot the position, velocity, and acceleration of the base vs time
sliderPosition = sys.myBodies{5}.myRTotal;
sliderVelocity = sys.myBodies{5}.myRDotTotal;
sliderAccel = sys.myBodies{5}.myRDDotTotal;
time = sys.myBodies{4}.myTimeTotal;

figure
subplot(3,1,1)
hold on
plot(time,sliderPosition(1,:))
plot(time,sliderPosition(2,:))
plot(time,sliderPosition(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Position (m)')
title('Slider')
hold off


subplot(3,1,2)
hold on
plot(time,sliderVelocity(1,:))
plot(time,sliderVelocity(2,:))
plot(time,sliderVelocity(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
title('Slider')
hold off


subplot(3,1,3)
hold on
plot(time,sliderAccel(1,:))
plot(time,sliderAccel(2,:))
plot(time,sliderAccel(3,:))
% axis([0 1 -3 5])
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Slider')
hold off


%% Plot the position, velocity, and acceleration of one side rod vs time
rodPosition = sys.myBodies{3}.myRTotal;
rodVelocity = sys.myBodies{3}.myRDotTotal;
rodAccel = sys.myBodies{3}.myRDDotTotal;
time = sys.myBodies{3}.myTimeTotal;

figure
subplot(3,1,1)
hold on
plot(time,rodPosition(1,:))
plot(time,rodPosition(2,:))
plot(time,rodPosition(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Position (m)')
title('Rod')
hold off


subplot(3,1,2)
hold on
plot(time,rodVelocity(1,:))
plot(time,rodVelocity(2,:))
plot(time,rodVelocity(3,:))
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
title('Rod')
hold off


subplot(3,1,3)
hold on
plot(time,rodAccel(1,:))
plot(time,rodAccel(2,:))
plot(time,rodAccel(3,:))
% axis([0 1 -3 5])
legend('x','y','z')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
title('Rod')
hold off








