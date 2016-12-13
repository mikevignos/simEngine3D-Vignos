function [ y, iterationCount ] = fourthOrderBDF( f, fJacobian, y0Vec, t0, tEnd, stepSize )
%backwardEuler.m
% Perform backward Euler for a giving intial value problem
%
% Function inputs:
% f : function
%   Function for first time derivative of y.
%
% gJacobian : function
%   Function for the jacobian of f.
%
% y0Vec : 4x1 double
%   Vector containing the values for y0, y1, y2, and y3. Since we are using
%   BDF 4th order, we need the first 4 terms of y provided.
%
% t0 : double
%   Start time for analysis.
%
% tEnd : double
%   Ending time for analysis.
%
% stepSize : double
%   Integration time step size.

% Set tolerance and max iterations for convergence of newton-raphson method
tolerance = 1e-8;
maxIter = 50;

% Initialize time vector and y
time = t0:stepSize:tEnd;
y = zeros(1,length(time));
iterationCount = zeros(1,length(time));
y(1:4) = y0Vec;

% Perform backward Euler method. Start from n = 5 since we are given
% y(1:4).
for iT = 5:length(time)
    tn = time(iT);
    
    % Initialize Newton-Raphson method
    yIter = y(iT - 1);
    ynMinus1 = y(iT - 1);
    ynMinus2 = y(iT - 2);
    ynMinus3 = y(iT - 3);
    ynMinus4 = y(iT - 4);
    iter = 1;
    while iter < maxIter
        % Compute Jacobian and g for this iteration
        gVal = yIter - 48/25*ynMinus1 + 36/25*ynMinus2 - 16/25*ynMinus3 + 3/25*ynMinus4 - 12/25*stepSize*f(tn,yIter);
        jacobian = 1 - 12/25*stepSize*fJacobian(tn,yIter);
        
        % Compute correction factor
        correction = jacobian\-gVal;
        
        % Apply correction factor
        yIter = yIter + correction;
        
        % Check if norm of residual is below the specified tolerance
        if norm(correction) < tolerance
            break;
        end
        iter = iter + 1;
    end
    iterationCount(iT) = iter;
    y(iT) = yIter;
end

end

