%% fourBarMechanismDynamicsAnalysis.m
% Taken from http://lim.ii.udc.es/mbsbenchmark/dist/A02/A02_specification.xml
% Used to test the combinaton of a revolute, universal, and spherical joint.  .
clear; close all; clc;

%% Create instance of multibody system
sys = multibodySystem();

%% Define parameters for the N 4-bar mechanism system
num4bar = 2; % Number of 4 bar mechanisms
L = 1; % m. Length of each link
mass = 1; % kg. Mass of each link.
v0 = 1; % m/s. Initial velocity of point B0

% Define inertial properties of each link
Jxx = mass*(L)^2/12;
Jyy = 0.0;
Jzz = mass*(L)^2/12;
JMatrix = zeros(3,3);
JMatrix(1,1) = Jxx;
JMatrix(2,2) = Jyy;
JMatrix(3,3) = Jzz;

%% Define the ground
mass1 = 1;
length1 = 0;
isGround1 = 1;
JMatrix1 = eye(3,3);
gravityDirection = '-y';
sys.addBody(1, 'ground', isGround1, mass1, length1, JMatrix1, gravityDirection);

%% Define all links in system
isGround = 0;
for iB = 1:(2*num4bar + 1)
    sys.addBody((iB+1), 'bar', isGround, mass, L, JMatrix, gravityDirection);
end

%% Define revolute joint between each horizontal link and its neighboring vertical links
% Hint: all positive integer bodies are the vertical links and all negative
% integer bodies are the horizontal links (excluding 1, which is the
% ground)
for iB = 3:2:(2*num4bar + 1)
    a1.body1 = iB - 1;
    a1.body2 = iB;
    a1.pointOnBody1 = [0 L/2 0]';
    a1.pointOnBody2 = [0 L/2 0]';
    a1.vector1OnBody1 = [1 0 0]';
    a1.vector2OnBody1 = [0 1 0]';
    a1.vectorOnBody2 = [0 0 1]';
    a1.constraintName = ['Rev joint b/w link ' num2str(iB - 1) ' and link ' num2str(iB)];
    sys.addJoint('revolute',a1);
    
    a2.body1 = iB + 1;
    a2.body2 = iB;
    a2.pointOnBody1 = [0 L/2 0]';
    a2.pointOnBody2 = [0 -L/2 0]';
    a2.vector1OnBody1 = [1 0 0]';
    a2.vector2OnBody1 = [0 1 0]';
    a2.vectorOnBody2 = [0 0 1]';
    a2.constraintName = ['Rev joint b/w link ' num2str(iB) ' and link ' num2str(iB+1)];
    sys.addJoint('revolute',a2);
end

%% Define revolute joints between each vertical link and the ground
% Hint: all positive integer bodies are the vertical links and all negative
% integer bodies are the horizontal links (excluding 1, which is the
% ground)
count = 1;
for iB = 2:2:(2*num4bar + 2)
    a.body1 = 1;
    a.body2 = iB;
    xPos = iB - (count + 1);
    a.pointOnBody1 = [xPos 0 0]';
    a.pointOnBody2 = [0 -L/2 0]';
    a.vector1OnBody1 = [1 0 0]';
    a.vector2OnBody1 = [0 1 0]';
    a.vectorOnBody2 = [0 0 1]';
    a.constraintName = ['Rev joint b/w ground and link ' num2str(iB)];
    sys.addJoint('revolute',a);
    count = count + 1;
end

%% Set inital position and orientation for all bodies
% if exist('NfourBarMechanismDynamicsAnalysis.mat')
%     load('NfourBarMechanismDynamicsAnalysis.mat')
% else
    rInitial = zeros(3,(2*num4bar + 2));
    
    pInitial =zeros(4,(2*num4bar + 2));
    pInitial(:,1) = [1 0 0 0]';
    pVert = [1 0 0 0]';
    Ahorz = [0 -1 0;
        1 0 0;
        0 0 1];
    pHorz = simEngine3DUtilities.A2p(Ahorz);
    
    count = 1;
    for iB = 2:(2*num4bar + 2)
        if (mod(iB,2) == 1) % Odd number
            xPos = (iB - 2)*L/2;
            yPos = L;
            zPos = 0;
            pTemp = pHorz(:);
            
        else % Even number
            
            xPos = iB - (count + 1);
            yPos = L/2;
            zPos = 0;
            pTemp = pVert(:);
            
            count = count + 1;
        end
        rInitial(:,iB) = [xPos yPos zPos]';
        pInitial(:,iB) = pTemp;
    end
    nBodies = (2*num4bar + 2);
    rInitial = reshape(rInitial,[3*nBodies,1]);
    pInitial = reshape(pInitial,[4*nBodies,1]);
    
    assemblyAnalysisFlag = 1;
    sys.setInitialPose( rInitial, pInitial, assemblyAnalysisFlag);
    
%     save('NfourBarMechanismDynamicsAnalysis.mat','sys')
% end

% Plot starting configuration of system
sys.plot(1);
    
%% Set inital velocities for all bodies
known = 2;
knownInitialRDot = [v0/2 0 0]';
knownInitialOmegaBar = [0; 0; -v0/L];
knownInitialPDot = simEngine3DUtilities.omegaBar2pDot(sys, known, knownInitialOmegaBar);

sys.computeAndSetInitialVelocities(known, knownInitialRDot, knownInitialPDot);

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
    save('NfourBarMechanismDynamicsAnalysis.mat','sys');
else
    load('NfourBarMechanismDynamicsAnalysis.mat')
end

disp(['Dynamics Analysis for NfourBarMechanismDynamicsAnalysis took ' num2str(analysisTime) ' seconds.'])

%% Plot the position of point B0 versus time
linkOrientation = sys.myBodies{2}.myPTotal;
linkPosition = sys.myBodies{2}.myRTotal;
time = 0:0.01:3.75;%sys.myBodies{2}.myTimeTotal;
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
