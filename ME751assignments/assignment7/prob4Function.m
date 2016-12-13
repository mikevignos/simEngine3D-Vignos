function [ g, J ] = prob4Function( xn, yn, xnMinus, ynMinus, h, alpha, beta )
% prob4Function.m

% Compute g(xn,yn) for the current guess of xn and yn
g(1,1) = xn*(1+h) + 4*h*xn*yn/(1 + xn^2) - xnMinus - h*alpha;
g(2,1) = -h*beta*xn + yn + h*beta*xn*yn/(1 + xn^2) - ynMinus;

% Compute jacobian for current guess of xn and yn
J(1,1) = 1 + h + 4*h*yn*(1 - xn^2)/(1 + xn^2)^2;
J(2,1) = -h*beta + h*beta*(1 - xn^2)/(1 + xn^2)^2;
J(1,2) = 4*h*xn/(1 + xn^2);
J(2,2) = 1 + beta*h*xn/(1 + xn^2);

end

