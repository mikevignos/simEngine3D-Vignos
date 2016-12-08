function bodyTrace(freeBodyNumber,state)
% simple comet animation of one body of system3D.
% INPUTS:
%   state : cell array of system states throughout time, this is generated
%           by kinematicsAnalysis, inverseDynamicsAnalysis, and dynamicsAnalysis functions

figure()
axis square

% allocate for body handles
x = zeros(length(state),1);
y = zeros(length(state),1);
z = zeros(length(state),1);

for i = 1:length(state)
    j = freeBodyNumber;
    x(i) = state{i}.r(3*j-2);
    y(i) = state{i}.r(3*j-1);
    z(i) = state{i}.r(3*j);
end
comet3(x,y,z);
end