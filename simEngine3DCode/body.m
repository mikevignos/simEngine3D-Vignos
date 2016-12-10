classdef body < handle
    %body.m Class that defines the attributes of a body in the system
    %   Detailed explanation goes here
    
    properties
        myBodyNumber; % Number of body
        myBodyType; % Currently only bars and ground are supported
        myIsGround; % Flag if the body is the ground
        myMass; % Mass of the body
        myMassMatrix; % Mass matrix for this body.
        myLength; % Length of the body. This will only be set if the body is a bar.
        myJMatrix; % Polar momnet of inertia matrix for this body.
        myJpMatrix; % 4*G'*J*G
        myG; % Current G matrix for the body
        myGDot; % Current Gdot matrix for this body.
        myGDDot; % Current Gddot matrix for this body.
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
        myRTotal; % Position of body across all time steps
        myRDotTotal; % Velocity of body across all time steps
        myRDDotTotal; % Acceleration of body across all time steps
        myPTotal; % Euler parameters of body across all time steps
        myPDotTotal; % First time derivate of euler parameters of body across all time steps
        myPDDotTotal; % Second time derivate of euler parameters of body across all time steps
        myForces; % Forces applied to this body.
        myConstraintForces; % Forces applied to this body due to each constraint in system.
        myConstraintTorques; % Torques applied to this body due to each constraint in system.
        myConstraintTorquesOmega; % Torques in r-omega formulation applied to this body due to each constraint in system.
        myConstraintForcesTotal; % Forces applied to this body due to each constraint in system throughout the entire simulation
        myConstraintTorquesTotal; % Torques applied to this body due to each constraint in system throughout the entire simulation.
        myConstraintTorquesOmegaTotal; % Torques in r-omega formulation applied to this body due to each constraint in system throughout the entire simulation.
        myTorques; % Torques applied to this body.
        myTorqueHat; % torqueHat vector for this specific body. This value 
                     % depends on the current orientation of the body, so 
                     % it is only computed when the user requests it.
        myTorqueHatPartialP; % Partial derivative of torqueHat vector w.r.t. p
        myTorqueHatPartialPDot; % Partial derivative of torqueHat vector w.r.t. pDot
        myJpDDotPartialP; % Partial derivate of [Jp*pDDot] w.r.t. p for this body.
    end
    
    properties (Dependent)
        myNumPoints; % Number of points defined on body
        myNumVectors; % Number of vectors defined on body
        myNumTimeSteps; % Total number of time steps
        myNumForces ; % Number of active forces being applied to this body.
        myNumTorques; % Number of active toruqes being appliced to this body.
        myTotalForce; % 3x1 vector containing the sum of all the forces applied to the body.
    end
    
    methods
        function obj = body(bodyNumber, bodyType, isGround, mass, bodyLength, JMatrix, gravityDirection)
            % Set parameters that are sent in
            obj.myBodyNumber = bodyNumber;
            obj.myBodyType = bodyType;
            obj.myIsGround = isGround;
            obj.myMass = mass;
            obj.myLength = bodyLength;
            
            % Store the user input Jmatrix of the body. It would be nice to
            % compute the Jmatrix for the user, but for most bodies J
            % depends on the choice of local reference frame.
            obj.myJMatrix = JMatrix;
            
            % Compute mass matrix and polar moment of inertia matrix for
            % this body.
            obj.myMassMatrix = mass*eye(3,3);
            
            % Add the force of gravity to this body. This currently assumes
            % that gravity acts in the -Z direction.
            if (isGround == 0)
                if strcmp(gravityDirection,'-z')
                    force = [0 0 -9.81*mass]';
                elseif strcmp(gravityDirection,'-y')
                    force = [0 -9.81*mass 0]';
                else
                    error('This gravity direction currently not implemented.');
                end
                sBar = [0 0 0]';
                obj.addForce(force, sBar, 'Force of Gravity');
            end
            
            % Set intial values
            obj.myRTotal = zeros(3,1); % Position of body across all time steps
            obj.myRDotTotal = zeros(3,1); % Velocity of body across all time steps
            obj.myRDDotTotal = zeros(3,1); % Acceleration of body across all time steps
            obj.myPTotal = zeros(4,1); % Euler parameters of body across all time steps
            obj.myPDotTotal = zeros(4,1); % First time derivate of euler parameters of body across all time steps
            obj.myPDDotTotal = zeros(4,1);
            
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
        function obj = convertConstraintTorques(obj)
            % Convert the constraint torques from the r-p formulation to
            % the r-omega formulation.
            
            torque = obj.myConstraintTorques;
            
            % Compute G for this body.
            obj.computeGmatrix();
            G = obj.myG;
            
            % Seed the new torque matrix
            [~,nConst] = size(torque);
            torqueOmega = zeros(3,nConst);
            
            for iC = 1:nConst
                torqueOmega(:,iC) = 0.5*G*torque(:,iC);
            end
            obj.myConstraintTorquesOmega = torqueOmega;
        end
        function obj = computeGmatrix(obj)
            % Compute the current G matrix for the body
            p = obj.myP;
            e0 = p(1);
            e = p(2:4);
            eTilde = simEngine3DUtilities.skewSym(e);
            G2 = -eTilde + e0*eye(3,3);
            G = [-e, G2];
            obj.myG = G;            
        end
        function obj = computeGDotMatrix(obj)
            % Compute the current G matrix for the body
            pDot = obj.myPDot;
            e0Dot = pDot(1);
            eDot = pDot(2:4);
            eDotTilde = simEngine3DUtilities.skewSym(eDot);
            G2 = -eDotTilde + e0Dot*eye(3,3);
            GDot = [-eDot, G2];
            obj.myGDot = GDot;            
        end
        function obj = computeGDDotMatrix(obj)
            % Compute the current GDDot matrix for the body
            pDDot = obj.myPDDot;
            e0DDot = pDDot(1);
            eDDot = pDDot(2:4);
            eDDotTilde = simEngine3DUtilities.skewSym(eDDot);
            G2 = -eDDotTilde + e0DDot*[1 0 0; 0 1 0; 0 0 1];
            GDDot = [-eDDot, G2];
            obj.myGDDot = GDDot;  
        end
        function obj = computeJpMatrix(obj)
            % Extract J matrix for body.
            J = obj.myJMatrix;
            
            % Compute G matrix
            obj.computeGmatrix();
            G = obj.myG;
            
            % Compute Jp matrix.
            Jp = 4*G'*J*G;
            obj.myJpMatrix = Jp;            
        end
        function obj = addForce(obj, force, sBar, forceName)
            % Add a force to this body. If this force also contributes a
            % torque, that needs to be computing during the dynamics
            % analysis phase.
            % 
            % Function inputs:
            % force : 3x1 double
            %   x,y, and z components of force defined in the global
            %   reference frame.
            %
            % sBar : 3x1 double
            %   Point of application of the force in the body reference
            %   frame. If the force is applied to the center of mass 
            %   sBar = [0 0 0]'.
            %
            % forceName : string
            %   Name of this force. Optional inpur.
            %   
            sBar = sBar(:);
            force = force(:);
            nForces = obj.myNumForces;
            bodyNum = obj.myBodyNumber;
            if nargin < 4
                forceName = ['Force ' num2str(nForces+1) ' on body ' num2str(bodyNum) ];
            end
            
            obj.myForces{nForces + 1}.force = force;
            obj.myForces{nForces + 1}.pointOfApplication = sBar;
            obj.myForces{nForces + 1}.name = forceName;
        end
        function obj = addTorque(obj, torque, torqueName)
            % Add a torque to this body.
            %
            % Function inputs:
            % torque : 3x1 double
            %   x, y, and z-components of torque.
            %
            % torqueName : string
            %   Name of this torque.
            
            torque = torque(:);
            nTorques = obj.myNumTorques;

            if nargin < 3
                torqueName = ['Torque ' num2str(nTorques+1) ' on body ' num2str(obj.myBodyNumber)];
            end
            
            obj.myTorques{nTorques + 1}.torque = torque;
            obj.myTorques{nTorques + 1}.name = torqueName;
        end
        function obj = computeTorqueHat(obj, sys, time)
            % Compute torqueHat vector for this body. This depends on the
            % current orientation of the body.
            
            % Compute the sum of all pure torques.
            torqueSum = [0 0 0]';
            for iT = 1:obj.myNumTorques
                torque = obj.myTorques{iT}.torque;
                torqueSum = torqueSum + torque;                
            end
            
            % Compute the sum of all torques due to forces. If the force is
            % applied at the center of mass, the variable torque (computed
            % below) will just be [0 0 0]'
            for iF = obj.myNumForces
                sBar = obj.myForces{iF}.pointOfApplication;
                force = obj.myForces{iF}.force;
                
                % Compute contribution of the force to the torque.
                    sBarTilde = simEngine3DUtilities.skewSym(sBar);
                    obj.computeA();
                    A = obj.myA;
                    torque = sBarTilde*A'*force;
                    torqueSum = torqueSum + torque;                
            end
            
            % Check if any TSDAs exist in the system. If they do, check to 
            % see if they are attached to this body. If they are, loop
            % through those and compute the torque due to the TSDA. Then,
            % add this into the total torque sum. Remember,
            % a torque is only applied to the body if the TSDA does not act
            % through the center of mass of the body.
            nTSDAs = sys.myNumTSDAs;
            if (nTSDAs > 0)
                for iT = 1:nTSDAs
                    TSDA = sys.myTSDAs{iT};
                    bodyI = TSDA.bodyI;
                    bodyJ = TSDA.bodyJ;
                    
                    % Only apply a torque if this body is one of the bodies
                    % the TSDA is attached to.
                    if (bodyI == obj.myBodyNumber) || (bodyJ == obj.myBodyNumber)
                        sBarIP = TSDA.sBarIP;
                        sBarJQ = TSDA.sBarJQ;
                        k = TSDA.stiffness;
                        L0 = TSDA.restingLength;
                        c = TSDA.dampingCoefficient;
                        h = TSDA.actuatorFunction;
                        
                        % Compute the force of the TSDA.
                        dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                        dijDot = simEngine3DUtilities.computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                        
                        lij = sqrt(dij'*dij);
                        eij = dij/lij;
                        lijDot = eij'*dijDot;
                        
                        force = k*(lij - L0) + c*lijDot + h(lij, lijDot, time);
                        
                        % Change direction of force depending on which body
                        % this is
                        if (bodyI == obj.myBodyNumber)
                            forceOnThisBody = force*eij;
                            sBar = sBarIP;
                        elseif (bodyJ == obj.myBodyNumber)
                            forceOnThisBody = -force*eij;
                            sBar = sBarJQ;
                        end
                        
                        % Compute the torque applied by this body.
                        obj.computeA();
                        A = obj.myA;
                        sBarTilde = simEngine3DUtilities.skewSym(sBar);
                        torque = sBarTilde*A'*forceOnThisBody;
                        
                        % Add this torque to the overall torque sum
                        torqueSum = torqueSum + torque;
                    end
                end
            end
            
            % Compute G and Gdot
            obj.computeGmatrix();
            obj.computeGDotMatrix();
            G = obj.myG;
            GDot = obj.myGDot;
            
            % Extract polar moment of inertia matrix and Euler params for
            % this body.
            p = obj.myP;
            Jbar = obj.myJMatrix;
            
            % Compute torqueHat for this body.
            torqueHat = 2*G'*torqueSum + 8*GDot'*Jbar*GDot*p;
            obj.myTorqueHat = torqueHat; 
        end
        function obj = computeJpDDotPartialP(obj)
            % Compute partial derivative of torqueHat w.r.t pDot for this
            % body.
            % Compute G, Gdot, Gddot
            obj.computeGmatrix();
            G = obj.myG;
            
            obj.computeGDotMatrix();
            GDot = obj.myGDot;
            
            obj.computeGDDotMatrix();
            GDDot = obj.myGDDot;
            
            % Extract polar moment of inertia matrix and p
            pDDot = obj.myPDDot;
            Jbar = obj.myJMatrix;
            
            % Compute JpDDotPartialP for this body
            a = Jbar*G*pDDot;
            aTilde = simEngine3DUtilities.skewSym(a);
            Ta = [0 -a';
                a -aTilde];
            JpDDotPartialP = 4*(Ta - G'*Jbar*GDDot);
            obj.myJpDDotPartialP = JpDDotPartialP;
        end
        function obj = computeTorqueHatPartialPDot(obj)
            % Compute partial derivative of torqueHat w.r.t pDot for this
            % body.
            % Compute G and Gdot
            obj.computeGmatrix();
            G = obj.myG;
            
            obj.computeGDotMatrix();
            GDot = obj.myGDot;
            
            % Extract polar moment of inertia matrix and p
            p = obj.myP;
            Jbar = obj.myJMatrix;
            
            % Compute partial derivative of torqueHat w.r.t pDot for this
            % body.
            a = Jbar*GDot*p;
            aTilde = simEngine3DUtilities.skewSym(a);
            Ta = [0 -a';
                a -aTilde];
            torqueHatPartialPDot = -8*GDot'*Jbar*G + 8*Ta;
            
            obj.myTorqueHatPartialPDot = torqueHatPartialPDot;
            
        end
        function obj = computeTorqueHatPartialP(obj)
            % Compute partial derivative of torqueHat w.r.t p for this
            % body.
            % Compute Gdot
            obj.computeGDotMatrix();
            GDot = obj.myGDot;
            
            % Extract polar moment of inertia matrix and Euler params for
            % this body.
            Jbar = obj.myJMatrix;
            
            % Compute torqueHat for this body.
            torqueHatPartialP = 8*GDot'*Jbar*GDot;
            
            obj.myTorqueHatPartialP = torqueHatPartialP;
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
            B = 2*[B1, B2];
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
        function myNumForces = get.myNumForces(obj)
            myNumForces = length(obj.myForces);
        end
        function myNumTorques = get.myNumTorques(obj)
            myNumTorques = length(obj.myTorques);
        end
        function myTotalForce = get.myTotalForce(obj)
            % Compute sum of all the forces acting on this body each time a
            % new force is added to the system. If no forces are added to
            % the system, myTotalForce = [0 0 0
            nForces = obj.myNumForces;
            forceTotal = zeros(3,1);
            for iF = 1:nForces
                forceTotal = forceTotal + obj.myForces{iF}.force;
            end
            myTotalForce = forceTotal;           
        end
    end
    
end
