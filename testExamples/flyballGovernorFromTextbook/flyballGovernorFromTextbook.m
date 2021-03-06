%% flyballGovernorFromTextbook.m
% Taken from ch. 12.6 section Edward J. Haug: Computer Aided Kinematics and Dynamics of
% Mechanical Systems (Allyn and Bacon, 1989)
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

% Add body2. The is the rotation shaft.
length2 = 0.2; % meters.
mass2 = 200;

% Define inertial properties of the rotation shaft.
Jxx = 25;
Jyy = 50;
Jzz = 25;
JMatrix2 = zeros(3,3);
JMatrix2(1,1) = Jxx;
JMatrix2(2,2) = Jyy;
JMatrix2(3,3) = Jzz;

isGround2 = 0;
sys.addBody(2, 'bar', isGround2, mass2, length2, JMatrix2, gravityDirection);

% Add body3 to the system. This is the right side ball.
length3 = 0.0; % meters.
mass3 = 1;

% Define inertial properties of the rotation shaft.
Jxx = 0.1;
Jyy = 0.1;
Jzz = 0.1;
JMatrix3 = zeros(3,3);
JMatrix3(1,1) = Jxx;
JMatrix3(2,2) = Jyy;
JMatrix3(3,3) = Jzz;

isGround3 = 0;
sys.addBody(3, 'point', isGround3, mass3, length3, JMatrix3, gravityDirection);

% Add body4 to the system. This is the right side ball in the picture.
length4 = 0.0; % meters.
mass4 = 1;

% Define inertial properties
Jxx = 0.1;
Jyy = 0.1;
Jzz = 0.1;
JMatrix4 = zeros(3,3);
JMatrix4(1,1) = Jxx;
JMatrix4(2,2) = Jyy;
JMatrix4(3,3) = Jzz;

isGround4 = 0;
sys.addBody(4, 'point', isGround4, mass4, length4, JMatrix4, gravityDirection);

% Add body5. This is the collar.
length5 = 0.1; % This length doesn't matter
mass5 = 1;

% Define inertial properties of the slider
Jxx = 0.15;
Jyy = 0.125;
Jzz = 0.15;
JMatrix5 = zeros(3,3);
JMatrix5(1,1) = Jxx;
JMatrix5(2,2) = Jyy;
JMatrix5(3,3) = Jzz;

isGround5 = 0;
sys.addBody(5, 'block', isGround5, mass5, length5, JMatrix5, gravityDirection);

%% Add TSDA to the system
s1.bodyI = 2;
s1.bodyJ = 5;
s1.sBarIP = [0 0 0]';
s1.sBarJQ = [0.0 0 0]';
s1.stiffness = 1000;
s1.restingLength = 0.15;
s1.dampingCoefficient = 30;
s1.actuatorFunction = @(lij, lijDot, t)0;
s1.name = 'Spring-Damper Element b/w body 2 and 5';
sys.addTSDA(s1);

%% Add a variable torques to apply to the axis.
bodyNumber = 2;
torqueFunction = @Ts;
torqueName = 'Ts';
sys.addVariableTorque(bodyNumber, torqueFunction, torqueName);

bodyNumber = 2;
torqueFunction = @Te;
torqueName = 'Te';
sys.addVariableTorque(bodyNumber, torqueFunction, torqueName);


%% Define important points and vectors in the system
sys.addPoint(2, [0 -0.2 0]', 'base');
sys.addPoint(2, [0 0 0]', 'top');
sys.addPoint(3, [-0.16 0 0]', 'base');
% sys.addPoint(3, [0 0 0.5]', 'top');
sys.addPoint(4, [0.16 0 0]', 'base');
% sys.addPoint(5, [0.05 0 0]','connectionForTSDA1');
% sys.addPoint(5, [-0.05 0 0]','connectionForTSDA2');

%% Define revolute joint between axis and ground
a1.body1 = 1;
a1.body2 = 2;
a1.pointOnBody1 = [0 0 0]';
a1.pointOnBody2 = [0 -0.2 0]';
a1.vector1OnBody1 = [1 0 0]';
a1.vector2OnBody1 = [0 0 1]';
a1.vectorOnBody2 = [0 1 0]';
a1.constraintName = 'Revolute Joint b/w Ground and Axis';

sys.addJoint('revolute',a1);

%% Define translational joint between axis and collar
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
a2.pointOnBody1 = [0 1 0]';
a2.pointOnBody2 = [0 0 0]';
a2.vector1OnBody1 = [1 0 0]';
a2.vector2OnBody1 = [0 0 1]';
a2.vector1OnBody2 = [0 1 0]';
a2.vector2OnBody2 = [0 0 1]';
a2.constraintName = 'Translational joint b/w axis and collar';

sys.addJoint('translational',a2);

%% Define revolute joint between axis and right ball
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
a3.pointOnBody1 = [0.0 0 0]';
a3.pointOnBody2 = [-0.16 0 0]';
a3.vector1OnBody1 = [1 0 0]';
a3.vector2OnBody1 = [0 1 0]';
a3.vectorOnBody2 = [0 0 1]';
a3.constraintName = 'Revolute Joint b/w right ball and Axis';

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
a4.pointOnBody1 = [0.0 0 0]';
a4.pointOnBody2 = [0.16 0 0]';
a4.vector1OnBody1 = [1 0 0]';
a4.vector2OnBody1 = [0 1 0]';
a4.vectorOnBody2 = [0 0 1]';
a4.constraintName = 'Revolute Joint b/w left ball and Axis';

sys.addJoint('revolute',a4);

%% Add distance constraint between right ball and collar
a5.bodyI = 3;
a5.bodyJ = 5;
a5.sBarIP = [-0.08 0 0]';
a5.sBarJQ = [0 0 0]';
a5.ft = @(t)0.10922^2;
a5.ftDot = @(t)0;
a5.ftDDot = @(t)0;
a5.constraintName = 'dist constraint b/w ball and collar';

isKinematic = 1;
sys.addBasicConstraint(isKinematic,'d',a5);

%% Add distance constraint between left ball and collar.
a6.bodyI = 4;
a6.bodyJ = 5;
a6.sBarIP = [0.08 0 0]';
a6.sBarJQ = [0 0 0]';
a6.ft = @(t)0.10922^2;
a6.ftDot = @(t)0;
a6.ftDDot = @(t)0;
a6.constraintName = 'dist constraint b/w ball and collar';

isKinematic = 1;
sys.addBasicConstraint(isKinematic,'d',a6);

%% Add driving constraint to model
% DP1 constraint between X axis of ground and z axis of shaft
% This is just for testing.
% a7.bodyJ = 1;
% a7.bodyI = 2;
% a7.aBarJ = [1 0 0]';
% a7.aBarI = [0 0 1]';
% a7.ft = @(t)cos(11.0174*t + pi/2 + 0.0001);
% a7.ftDot = @(t)(-11.0174*sin(pi/2 + 11.0174*t + 0.0001));
% a7.ftDDot = @(t)(-(11.0174)^2*cos(pi/2 + 11.0174*t + 0.0001));
% a7.constraintName = 'DP1 driving constraint';
% isKinematic = 0;
% sys.addBasicConstraint(isKinematic,'dp1',a7);

%% Set initial conditions of each body
% if exist('flyballGovernorSetup.mat')
%     load('flyballGovernorSetup.mat')
% else
% alpha = 30;
% r1Initial = zeros(3,1); % Ground
% r2Initial = [0; 0.2; 0]; % Axis
% r3Initial = [0.16*sind(45); 0.2 - 0.16*cosd(45); 0]; % Right rod
% r4Initial = [-0.16*sind(45); 0.2 - 0.16*cosd(45); 0]; % Left rod
% r5Initial = [0; 0.05; 0]; % Base
r1Initial = zeros(3,1); % Ground
r2Initial = [0; 0.2; 0]; % Axis
r3Initial = [0.11314; 0.08686; 0]; % Right rod
r4Initial = [-0.11314; 0.08686; 0]; % Left rod
r5Initial = [0; 0.05; 0]; % Base
rInitial = [r1Initial; r2Initial; r3Initial; r4Initial; r5Initial];

% Initial orientation
% Ground
p1 = [1.0 0.0 0.0 0.0]';

% Axis
p2 = [1 0 0 0]';

% Right rod
% s45 = sind(45);
% c45 = cosd(45);
% A3 = [c45 c45 0;
%     -s45 s45 0;
%     0 0 1];
% p3 = simEngine3DUtilities.A2p(A3);
p3 = [-0.9239 0 0 0.3827]';

% Left rod
% % A4 = [c45 -c45 0;
% %     s45 s45 0;
% %     0 0 1];
% % p4 = simEngine3DUtilities.A2p(A4);
p4 = [0.9239 0 0 0.3827]';

% Base
p5 = [1 0 0 0]';

pInitial = [p1; p2; p3; p4; p5];

t = 0;
assemblyAnalysisFlag = 1;
sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);

% Initial velocities. Initial velocity known for crank.
known = 2;
knownInitialRDot = [0; 0; 0];
knownInitialOmegaBar = [0; 11.0174; 0];
knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);
sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);

% Use the next command if you are prescribing a driving constraint.



% sys.computeAndSetInitialVelocities([], [], []);


%     save('flyballGovernorSetup.mat','sys');
% end

%% Plot starting configuration of system
sys.plot(1);
view([2])

%% Perform analysis
if 1
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-3;
    order = 2;
    displayFlag = 1;
    velocityConstraintFlag = 0;
    method = 'quasiNewton';
    tic;
    sys.dynamicsAnalysis(timeStart, timeEnd,timeStep, order, method, displayFlag, velocityConstraintFlag);
    analysisTime = toc;
    save('flyballGovernorFromText_k1000_C12500.mat','sys');
    disp(['Analysis for flyballGovernorFromText_k1000_C12500 took ' num2str(analysisTime) ' seconds.'])
else
    load('flyballGovernorFromText_k1000_C12500.mat')
end



%% Animate results
plot.animateSystem(sys,[0 90])

%% Plot the position, velocity, and acceleration of the base vs time
sliderPosition = sys.myBodies{5}.myRTotal;
sliderVelocity = sys.myBodies{5}.myRDotTotal;
sliderAccel = sys.myBodies{5}.myRDDotTotal;
time = sys.myBodies{5}.myTimeTotal;

% Plot zoomed in version of slider position
sliderPositionCM = sliderPosition*100;
figure
hold on
plot(time,sliderPositionCM(2,:),'LineWidth',2)
legend('x','y','z')
axis([timeStart timeEnd 4.6 5.1])
xlabel('Time (sec)','FontSize',16)
ylabel('Position (cm)','FontSize',16)
title('Slider Y-Position','FontSize',16)
set(gca, 'FontSize',12);
hold off

% figure
% subplot(3,1,1)
% hold on
% plot(time,sliderPosition(1,:))
% plot(time,sliderPosition(2,:))
% plot(time,sliderPosition(3,:))
% legend('x','y','z')
% xlabel('Time (sec)')
% ylabel('Position (m)')
% title('Slider')
% hold off
%
%
% subplot(3,1,2)
% hold on
% plot(time,sliderVelocity(1,:))
% plot(time,sliderVelocity(2,:))
% plot(time,sliderVelocity(3,:))
% legend('x','y','z')
% xlabel('Time (sec)')
% ylabel('Velocity (m/s)')
% title('Slider')
% hold off
%
%
% subplot(3,1,3)
% hold on
% plot(time,sliderAccel(1,:))
% plot(time,sliderAccel(2,:))
% plot(time,sliderAccel(3,:))
% % axis([0 1 -3 5])
% legend('x','y','z')
% xlabel('Time (sec)')
% ylabel('Acceleration (m/s^2)')
% title('Slider')
% hold off


%% Plot the angular velocity of the shaft versus time
pForShaftC7500 = sys.myBodies{2}.myPTotal;
pDotForShaft = sys.myBodies{2}.myPDotTotal;
time = sys.myBodies{2}.myTimeTotal;

% Convert pDot to omegaBar
omegaBar = zeros(3,length(time));
for iT = 1:length(time)
    % Compute the current G matrix for the body
    p = pForShaftC7500(:,iT);
    e0 = p(1);
    e = p(2:4);
    eTilde = simEngine3DUtilities.skewSym(e);
    G2 = -eTilde + e0*eye(3,3);
    G = [-e, G2];
    
    % Compute omageBar
    pDot = pDotForShaft(:,iT);
    omegaBar(:,iT) = 2*G*pDot;
end

% Plot omegaBar versus time
figure
hold on
plot(time,omegaBar(1,:));
plot(time,omegaBar(2,:));
plot(time,omegaBar(3,:));
axis([timeStart timeEnd 10.1 11.1])
xlabel('Time (sec)')
ylabel('\omega (rad/sec)');
legend('x','y','z')
title('Angular velocity of shaft')

%% Plotting of 3 different simulations
if 1
    sliderPositionC7500 = sysC7500.myBodies{5}.myRTotal*100;
    timeC7500 = sysC7500.myBodies{5}.myTimeTotal;
    
    sliderPositionC12500 = sysC12500.myBodies{5}.myRTotal*100;
    timeC12500 = sysC12500.myBodies{5}.myTimeTotal;
    
    sliderPositionC17500 = sysC17500.myBodies{5}.myRTotal*100;
    timeC17500 = sysC17500.myBodies{5}.myTimeTotal;
    
    % Plot zoomed in version of slider position
    figure
    hold on
    plot(timeC7500,sliderPositionC7500(2,:),'LineWidth',2)
    plot(timeC12500,sliderPositionC12500(2,:),'LineWidth',2)
    plot(timeC12500,sliderPositionC17500(2,:),'LineWidth',2)
    legend('C = 7500','C = 12500','C = 17500')
    axis([timeStart timeEnd 4.6 5.1])
    xlabel('Time (sec)','FontSize',16)
    ylabel('Position (cm)','FontSize',16)
    title('Slider Y-Position','FontSize',16)
    set(gca, 'FontSize',12);
    hold off
    
    % Plot the angular velocity of the shaft versus time
    pForShaftC7500 = sysC7500.myBodies{2}.myPTotal;
    pDotForShaftC7500 = sysC7500.myBodies{2}.myPDotTotal;
    timeC7500 = sysC7500.myBodies{2}.myTimeTotal;
    
    pForShaftC12500 = sysC12500.myBodies{2}.myPTotal;
    pDotForShaftC12500 = sysC12500.myBodies{2}.myPDotTotal;
    timeC12500 = sysC12500.myBodies{2}.myTimeTotal;
    
    pForShaftC17500 = sysC17500.myBodies{2}.myPTotal;
    pDotForShaftC17500 = sysC17500.myBodies{2}.myPDotTotal;
    timeC17500 = sysC17500.myBodies{2}.myTimeTotal;
    
    % Convert pDot to omegaBar
    omegaBarC7500 = zeros(3,length(timeC7500));
    for iT = 1:length(timeC7500)
        % Compute the current G matrix for the body
        p = pForShaftC7500(:,iT);
        e0 = p(1);
        e = p(2:4);
        eTilde = simEngine3DUtilities.skewSym(e);
        G2 = -eTilde + e0*eye(3,3);
        G = [-e, G2];
        
        % Compute omageBar
        pDot = pDotForShaftC7500(:,iT);
        omegaBarC7500(:,iT) = 2*G*pDot;
    end
    
    omegaBarC12500 = zeros(3,length(timeC12500));
    for iT = 1:length(timeC12500)
        % Compute the current G matrix for the body
        p = pForShaftC12500(:,iT);
        e0 = p(1);
        e = p(2:4);
        eTilde = simEngine3DUtilities.skewSym(e);
        G2 = -eTilde + e0*eye(3,3);
        G = [-e, G2];
        
        % Compute omageBar
        pDot = pDotForShaftC12500(:,iT);
        omegaBarC12500(:,iT) = 2*G*pDot;
    end
    
    omegaBarC17500 = zeros(3,length(timeC17500));
    for iT = 1:length(timeC17500)
        % Compute the current G matrix for the body
        p = pForShaftC17500(:,iT);
        e0 = p(1);
        e = p(2:4);
        eTilde = simEngine3DUtilities.skewSym(e);
        G2 = -eTilde + e0*eye(3,3);
        G = [-e, G2];
        
        % Compute omageBar
        pDot = pDotForShaftC17500(:,iT);
        omegaBarC17500(:,iT) = 2*G*pDot;
    end
    
    
    % Plot omegaBar versus time
    figure
    hold on
    plot(timeC7500,omegaBarC7500(2,:),'LineWidth',2);
    plot(timeC12500,omegaBarC12500(2,:),'LineWidth',2);
    plot(timeC17500,omegaBarC17500(2,:),'LineWidth',2);
    axis([timeStart timeEnd 10.1 11.1])
    xlabel('Time (sec)','FontSize',16)
    ylabel('\omega (rad/sec)','FontSize',16);
    legend('C = 7500','C = 12500','C = 17500')
    set(gca, 'FontSize',12);
    title('Angular Velocity of Shaft','FontSize',16)
    
end