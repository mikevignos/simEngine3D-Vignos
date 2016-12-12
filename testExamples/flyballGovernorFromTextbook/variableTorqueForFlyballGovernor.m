function [ torque ] = variableTorqueForFlyballGovernor( sys, time )
%Function that will be used to compute the torque to apply to the flyball
%governor.
C = 7500;

% Extract the current y-position of the collar
yLoc = sys.myBodies{5}.myR(2);
L0 = 0.15;
currentL = 0.2 - yLoc;
deltaL = currentL - L0;
torqueY = C*deltaL;
torque = [0; torqueY; 0];
end

