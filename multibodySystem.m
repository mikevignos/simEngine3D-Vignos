classdef multibodySystem < handle
    %SYSTEM.m
    %   Defines a multi-body system composed of rigid bodies, constraints,
    %   and joints.
    
    properties
        myBodies; % Structure containing all of the bodies in the system
        myConstraints; % Structure containing all of the constraints in the system
    end
    
    methods
        function obj = multibodySystem()
        end
        function obj = addBody(obj,bodyNumber, bodyType, isGround, mass, bodyLength)
            % mass : double
            %   Mass of body in kg
            %
            % bodyLength : double
            %   Length of body in meters
            
            % Add a body to the system
            newBody = body(bodyNumber, bodyType, isGround, mass, bodyLength);
            obj.myBodies{bodyNumber} = newBody;
        end
        function obj = addPoint(obj, bodyNumber, sBar, pointName)
            % Adds a point to this body. This function is helpful for
            % storing points that will be used in constraints.
            % 
            % Function inputs:
            % bodyNumber: double
            %   Number of body in this system.
            %
            % sBar : 3x1 vector
            %   Location of point in body reference frame
            %   
            % pointName : string
            %   Name for this point. This is not a necessary input
            %
            
            obj.myBodies{bodyNumber} = obj.myBodies{bodyNumber}.addPoint(sBar, pointName);
        end
        function obj = addVector(obj, bodyNumber, aBar, vectorName)
            % Adds a point to this body. This function is helpful for
            % storing points that will be used in constraints.
            % 
            % Function inputs:
            % bodyNumber: double
            %   Number of body in this system.
            %
            % aBar : 3x1 vector
            %   Vector in body reference frame
            %   
            % pointName : string
            %   Name for this point. This is not a necessary input
            %
            
            obj.myBodies{bodyNumber}.addVector(aBar, vectorName);
        end
        function obj = updateSystemState(obj, rMatrix, rDotMatrix, pMatrix, pDotMatrix, time)
            % Update position, orientation, and time derivatives of each
            % for each body at a specific time step. 
            %
            % Function inputs:
            % rMatrix : 3 x N matrix, where N is the number of bodies
            %   Matrix containing the current 3D positon of each body in
            %   the global reference frame
            %
            % rDotMatrix : 3 x N matrix, where N is the number of bodies
            %   Matrix containing the current time derivative of the 3D 
            %   positon of each body in the global reference frame
            %
            % pMatrix : 4 x N matrix, where N is the number of bodies
            %   Matrix containing the current Euler parameters of each body
            %
            % pDotMatrix : 4 x N matrix, where N is the number of bodies
            %   Matrix containing the current time derivative of the 
            %   Euler parameters of each body
            
            nBodies = length(obj.myBodies);
            for iB = 1:nBodies
                p = pMatrix(:,iB);
                pDot = pDotMatrix(:,iB);
                r = rMatrix(:,iB);
                rDot = rDotMatrix(:,iB);
                obj.myBodies{iB}.updateBody(p, pDot, r, rDot, time);
            end
        end
        function obj = addDP1constraint(obj,bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName)
            % Add a DP1 constraint between two bodies in the system
            newConstraint = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName);
            
            % Current number of constraints
            nConst = length(obj.myConstraints);
            obj.myConstraints{nConst + 1} = newConstraint;            
        end
        function obj = addCDconstraint(obj, bodyI, bodyJ, coordVec, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName)
            % Add a CD constraint between two bodies in the system
            newConstraint = CDconstraint(bodyI, bodyJ, coordVec, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName);
            
            % Current number of constraints
            nConst = length(obj.myConstraints);
            obj.myConstraints{nConst + 1} = newConstraint;    
        end
        function obj = computeConstraintProperties(obj, constraintNumber, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
            % Determine type of constraint
            constraintType = obj.myConstraints{constraintNumber}.myConstraintType;
            
            % Call the correct constraint class based off constraint type.
            switch constraintType
                case 'DP1'
                    obj.myConstraints{constraintNumber}.computeDP1constraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
                    
                case 'CD'
                    obj.myConstraints{constraintNumber}.computeCDconstraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
            end
        end
        
    end
    
end