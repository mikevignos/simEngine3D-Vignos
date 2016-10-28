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
a6.aBarI = [0 1 0]';
a6.ft = @(t)cos((pi*cos(2*t))/4 + pi/2);
a6.ftDot = @(t)((pi*sin(2*t)*sin((pi*cos(2*t))/4 + pi/2))/2);
a6.ftDDot = @(t)(pi*cos(2*t)*sin((pi*cos(2*t))/4 + pi/2) - (pi^2*sin(2*t)^2*cos((pi*cos(2*t))/4 + pi/2))/4);
a6.constraintName = 'DP1 driving constraint';
isKinematic = 0;
sys.addBasicConstraint(isKinematic,'dp1',a6);

%% Perform kinematics analysis
if 0
    timeStart = 0;
    timeEnd = 10;
    timeStep = 10^-3;
    sys.kinematicsAnalysis(timeStart, timeEnd, timeStep, 0);
    % Save the multibody system since it contains the results
    save('multibodySystem_A6P3.mat','sys');
    
else
    load('multibodySystem_A6P3.mat')
end


%% Create plots for origin of body 2 ref frame
time = sys.myBodies{2}.myTimeTotal;
rOprime = sys.myBodies{2}.myRTotal;
rDotOprime = sys.myBodies{2}.myRDotTotal;
rDDotOprime = sys.myBodies{2}.myRDDotTotal;
pOprime = sys.myBodies{2}.myPTotal;
pDotOprime = sys.myBodies{2}.myPDotTotal;
pDDotOprime = sys.myBodies{2}.myPDDotTotal;

figure
subplot(3,1,1)
hold on
plot(time,rOprime(1,:))
plot(time,rOprime(2,:))
plot(time,rOprime(3,:))
title('Position of point O-prime')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rDotOprime(1,:))
plot(time,rDotOprime(2,:))
plot(time,rDotOprime(3,:))
title('Velocity of point O-prime')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rDDotOprime(1,:))
plot(time,rDDotOprime(2,:))
plot(time,rDDotOprime(3,:))
title('Acceleration of point O-prime')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off

figure
hold on
plot(time,rDDotOprime(1,:))
plot(time,rDDotOprime(2,:))
plot(time,rDDotOprime(3,:))
axis([0 10 -6 6])
title('Acceleration of point O-prime -- Zoomed In')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off

%% Compute position, velocity, and acceleration of point Q
sBar = [-2 0 0]';
rQ = zeros(3,length(time));
rQdot = zeros(3,length(time));
rQddot = zeros(3,length(time));
for i = 1:length(time)
    p = pOprime(:,i);
    pDot = pDotOprime(:,i);
    pDDot = pDDotOprime(:,i);
    
    A = simEngine3DUtilities.p2A(p);
    B = simEngine3DUtilities.computeBmatrix(p,sBar);
    Bdot = simEngine3DUtilities.computeBmatrix(pDot,sBar);
    rQ(:,i) = rOprime(:,i) + A*sBar;
    rQdot(:,i) = rDotOprime(:,i) + B*pDot;
    rQddot(:,i) = rDDotOprime(:,i) + B*pDDot + Bdot*pDot;
end


%% Create plots for point Q
figure
subplot(3,1,1)
hold on
plot(time,rQ(1,:))
plot(time,rQ(2,:))
plot(time,rQ(3,:))
title('Position of point Q')
xlabel('Time (sec)')
ylabel('Position (m)')
legend('X','Y','Z')
hold off

subplot(3,1,2)
hold on
plot(time,rQdot(1,:))
plot(time,rQdot(2,:))
plot(time,rQdot(3,:))
title('Velocity of point Q')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
legend('X','Y','Z')
hold off

subplot(3,1,3)
hold on
plot(time,rQddot(1,:))
plot(time,rQddot(2,:))
plot(time,rQddot(3,:))
title('Acceleration of point Q')
xlabel('Time (sec)')
ylabel('Acceleration (m/s^2)')
legend('X','Y','Z')
hold off
