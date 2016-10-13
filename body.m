classdef body < handle
    %body.m Class that defines the attributes of a body in the system
    %   Detailed explanation goes here
    
    properties
        myBodyNumber; % Number of body
        myBodyType; % Currently only bars and ground are supported
        myIsGround; % Flag if the body is the ground
        myMass; % Mass of the body
        myLength; % Length of the body. This will only be set if the body is a bar.
        myJ; % Polar moment of inertia for body
        myPoints; % Structure containing the important points on the body.
        myVectors; % Structure containing important vectors on the body
        myP = [1; 0; 0; 0]; % Euler parameters of local reference frame for given body
        myPDot = [0; 0; 0; 0]; % Time derivative of Euler paramters
        myR = [0; 0; 0]; % Location of body in global reference frame
        myRDot = [0; 0; 0]; % Time derivation of body location (i.e. components of body velocity)
        myTime; % Current time step at which p, pDot, r, and rDot were last set.
        myA; % Orientation matrix for current euler paramters
        myB; % Current B matrix for body
        myBDot; % Time derivative of B matrix
    end
    
    methods
        function obj = body(bodyNumber, bodyType, isGround, mass, length)
            % Set parameters that are sent in
            obj.myBodyNumber = bodyNumber;
            obj.myBodyType = bodyType;
            obj.myIsGround = isGround;
            obj.myMass = mass;
            obj.myLength = length;
            
            % NEED TO WRITE. Compute polar moment of inertia for body
            
            
            
        end
        function obj = updateBody(obj, p, pDot, r, rDot, time)
            % Update the current orientation and position of the given body
            obj.myP = p;
            obj.myR = r;          
            obj.myPDot = pDot;
            obj.myRDot = rDot;
            obj.myTime = time;
        end
        function obj = addPoint(obj, sBar, pointName)
            % Adds a point to this body. This function is helpful for
            % storing points that will be used in constraints.
            % 
            % Function inputs:
            % sBar : 3x1 vector
            %   Location of point in body reference frame
            %   
            % pointName : string
            %   Name for this point. This is not a necessary input
            %
            
            % Determine how many points have already been created
            nPoints = length(obj.myPoints);
            
            % Set name of points
            if nargin < 3
                pointName = ['Point' num2str(nPoints + 1)];
            end

            % Make sure sBar is a column vector
            sBar = sBar(:);
            
            % Set properties of point
            obj.myPoints{nPoints + 1}.name = pointName;
            obj.myPoints{nPoints + 1}.sBar = sBar;           
        end
        function obj = addVector(obj, aBar, vectorName)
            % Adds a vector to this body. This function is helpful for
            % storing vectors that will be used in constraints.
            % 
            % Function inputs:
            % aBar : 3x1 vector
            %   Vector in body reference frame. If this is not a unit
            %   vector, the function will convert it to one.
            %   
            % pointName : string
            %   Name for this point. This is not a necessary input
            %
            
            % Determine how many vectors have already been created
            nVecs = length(obj.myVectors);
            
            % Set name of points
            if nargin < 3
                vectorName = ['Vector' num2str(nVecs + 1)];
            end

            % Make sure aBar is a column vector and a unit vector
            aBar = aBar(:);
            
            aBar = aBar/norm(aBar);
            
            % Set properties of vector
            obj.myVectors{nVecs + 1}.name = vectorName;
            obj.myVectors{nVecs + 1}.aBar = aBar;           
        end
        function obj = computeA(obj)
            % Compute the rotation matrix for this body in the current
            % orientation
            p = obj.myP;
            e0 = p(1);
            e1 = p(2); 
            e2 = p(3);
            e3 = p(4);
            
            A(1,1) = e0^2 + e1^2 - 0.5;
            A(1,2) = e1*e2 - e0*e3;
            A(1,3) = e1*e3 + e0*e2;
            A(2,1) = e1*e2 + e0*e3;
            A(2,2) = e0^2 + e2^2 - 0.5;
            A(2,3) = e2*e3 - e0*e1;
            A(3,1) = e1*e3 - e0*e2;
            A(3,2) = e2*e3 + e0*e1;
            A(3,3) = e0^2 + e3^2 - 0.5;
            Afinal = 2*A;
            
            obj.myA = Afinal;
        end
        function obj = computeB(obj,aBar)
            p = obj.myP;
            e0 = p(1);
            e = p(2:4);
            eTilde = skewSym(e);
            aBarTilde = skewSym(aBar);
            
            % Create B matrix
            B1 = (e0*eye(3,3) + eTilde)*aBar;
            B2 = e*aBar' - (e0*eye(3,3) + eTilde)*aBarTilde;
            B = 2*[B1 B2];
            obj.myB = B;
        end
        
        function obj = computeBDot(obj,aBar)
            pDot = obj.myPDot;
            e0Dot = pDot(1);
            eDot = pDot(2:4);
            
            eDotTilde = skewSym(eDot);
            aBarTilde = skewSym(aBar);
            
            % Compute BDot matrix
            B1 = (e0Dot*eye(3,3) + eDotTilde)*aBar;
            B2 = eDot*aBar' - (e0Dot*eye(3,3) + eDotTilde)*aBarTilde;
            BDot = 2*[B1 B2];
            obj.myBDot = BDot;      
        end       
    end
    
end
