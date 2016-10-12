classdef system
    %SYSTEM.m
    %   Defines a multi-body system composed of rigid bodies, constraints,
    %   and joints.
    
    properties
        myBodies; % Structure containing all of the bodies in the system
        myConstraints; % Structure containing all of the constraints in the system
    end
    
    methods
        function obj = system()
        end
        function obj = addBody(obj,bodyNumber, bodyType, mass, length)
            % Add a body to the system
            newBody = body(bodyNumber, bodyType, mass, length);
            obj.myBodies{bodyNumber} = newBody;
        end
        function obj = addDP1constraint(obj,bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName)
            % Add a DP1 constraint between two bodies in the system
            newConstraint = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName);
            
            % Current number of constraints
            nConst = length(obj.myConstraints);
            obj.myConstraints{nConst + 1} = newConstraint;            
        end
        function obj = addCDconstraint(obj)
        end
    end
    
end

