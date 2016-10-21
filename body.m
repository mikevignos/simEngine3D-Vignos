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
        myP; % Euler parameters of local reference frame for given body
        myPDot; % Time derivative of Euler paramters
        myPDDot; % Second time derivative of Euler parameters
        myR; % Location of body in global reference frame
        myRDot; % Time derivation of body location (i.e. components of body velocity)
        myRDDot; % Second time derivative of body location (i.e. components of body acceleration)
        myTime; % Current time step at which p, pDot, r, and rDot were last set.
        myA; % Orientation matrix for current euler paramters
        myB; % Current B matrix for body
        myBDot; % Time derivative of B matrix
        myTimeTotal; % Vector containing all time steps traversed in analysis
        myRTotal = zeros(3,1); % Position of body across all time steps
        myRDotTotal = zeros(3,1); % Velocity of body across all time steps
        myRDDotTotal = zeros(3,1); % Acceleration of body across all time steps
        myPTotal = zeros(4,1); % Euler parameters of body across all time steps
        myPDotTotal = zeros(4,1); % First time derivate of euler parameters of body across all time steps
        myPDDotTotal = zeros(4,1); % Second time derivate of euler parameters of body across all time steps
    end
    
    properties (Dependent)
        myNumPoints; % Number of points defined on body
        myNumVectors; % Number of vectors defined on body
        myNumTimeSteps; % Total number of time steps 
    end
    
    methods
        function obj = body(bodyNumber, bodyType, isGround, mass, bodyLength)
            % Set parameters that are sent in
            obj.myBodyNumber = bodyNumber;
            obj.myBodyType = bodyType;
            obj.myIsGround = isGround;
            obj.myMass = mass;
            obj.myLength = bodyLength;
            
            % Compute polar moment of inertia for body
            if strcmp(bodyType,'bar')
                J = (mass*bodyLength^2)/12;
            else
                J = 0;
            end
            obj.myJ = J;
            
            
        end
        
        function obj = updateBody(obj, p, pDot, pDDot, r, rDot, rDDot, time)
            % Update the current orientation and position of the given body
            obj.myP = p;
            obj.myR = r;          
            obj.myPDot = pDot;
            obj.myRDot = rDot;
            obj.myPDDot = pDDot;
            obj.myRDDot = rDDot;
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
            % If the make is the ground, then make sure the orientation
            % matrix never changes
            isGround = obj.myIsGround;
            if (isGround == 1)
                Afinal = eye(3,3);
            else
                % Compute the orientation matrix for this body in the current
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
            end
            
            
            obj.myA = Afinal;
        end
        function obj = computeB(obj,aBar)
            p = obj.myP;
            e0 = p(1);
            e = p(2:4);
            eTilde = simEngine3DUtilities.skewSym(e);
            aBarTilde = simEngine3DUtilities.skewSym(aBar);
            
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
            
            eDotTilde = simEngine3DUtilities.skewSym(eDot);
            aBarTilde = simEngine3DUtilities.skewSym(aBar);
            
            % Compute BDot matrix
            B1 = (e0Dot*eye(3,3) + eDotTilde)*aBar;
            B2 = eDot*aBar' - (e0Dot*eye(3,3) + eDotTilde)*aBarTilde;
            BDot = 2*[B1 B2];
            obj.myBDot = BDot;      
        end       
    end
    
    % Methods for dependent properties
    methods
        function myNumPoints = get.myNumPoints(obj)
            myNumPoints = length(obj.myPoints);
        end
        function myNumVectors = get.myNumVectors(obj)
            myNumVectors = length(obj.myVectors);
        end
        function myNumTimeSteps = get.myNumTimeSteps(obj)
            myNumTimeSteps = length(obj.myTimeTotal);
        end
    end
    
end
