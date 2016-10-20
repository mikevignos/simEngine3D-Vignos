classdef intermediateConstraint < handle
    % intermediateConstraint.m
    % Implements an intermediate geometric constraint. These constraints
    % are typically made up of a few basic constraints (DP1, DP2, CD, or D)
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType; % String. Type of constraint.
        myIsKinematic; % Flag for if the constraint is a kinematic (1) or driving (0) constraint
        myAttributes; % Attributes of the desired constraint. 
                      % The attributes will vary depending on which type of
                      % constraint is defined.
                      %%%%%%
%         myBodyJ; % Second body in constraint
%         mysBarIP; % Point on body I being constrained in body I reference frame
%         mysBarJQ; % Point on body J being constrained in body J reference frame
%         myFt; % Function representing the square of the distance between points mysBarIP and mysBarJQ
%         myFtDot; % Function representing first time derivative of f(t)
%         myFtDDot; % Function representing second time derivative of f(t)
        myTime; % Value of current time step
        
        % Properties of constraint that can be computed
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates (i.e. r)
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates (i.e. p)
    end
    
    methods
    end
    
end

