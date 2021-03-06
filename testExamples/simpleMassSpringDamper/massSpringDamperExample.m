clear; close all; clc;
tspan=[0 4];
y0=[0.02;0];
[t,y]=ode45('unforced1',tspan,y0);
plot(t,y(:,1));
grid on
xlabel('time')
ylabel('Displacement')
title('Mass-Spring-Damper Solved Using ODE45 in MATLAB')
hold on; 