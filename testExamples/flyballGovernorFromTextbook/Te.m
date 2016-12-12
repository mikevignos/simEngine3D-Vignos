function [ torque ] = Te( sys, time )
%Function that will be used to compute the torque to apply to the flyball
%governor.

% Extract the current y-position of the collar
if (time < 1)
    torqueY = 0;
elseif (time >=1) && (time < 2)
    torqueY = -25*time + 25;
else
    torqueY = -25;
end

torque = [0; torqueY; 0];

end

