%% ME 751 Homework 7 -- Vignos
clear; close all; clc;

%% Problem 3
% Run forward Euler for a second of simulation time. Plot the errors in the
% solution from the analytical solution at each timestep.
yDot = @(t,lambda)lambda*y;

% lambda = -10. This should go unstable at h = 0.2


% lambda = -100. This should go unstable at h = 0.02

% lambda = -1000. This should go unstable at h = 0.002