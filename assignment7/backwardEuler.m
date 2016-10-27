function [ y ] = backwardEuler( g, gJacobian, y0, t0, tEnd, stepSize )
%backwardEuler.m
% Perform backward Euler for a giving intial value problem
%
% Function inputs:
% g : function
%   Function for the newton iteration of the backward euler method. This
%   function is dependent on y(n), y(n-1), and time.
%
% gJacobian : function
%   Function representing the jacobian of g.
%
% y0 : double
%   Initial value for y.
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
y(1) = y0;

% Perform backward Euler method
for iT = 2:length(time)
    t = time(iT);
    
    % Initialize Newton-Raphson method
    yIter = y(iT - 1);
    ynMinus = y(iT - 1);
    iter = 1;
    while iter < maxIter
        % Compute Jacobian and g for this iteration
        gVal = g(yIter, ynMinus, stepSize, t);
        jacobian = gJacobian(yIter, stepSize,t);
        
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
    y(iT) = yIter;
end

end

