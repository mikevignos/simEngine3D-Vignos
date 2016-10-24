classdef multibodySystem < handle
    %SYSTEM.m
    %   Defines a multi-body system composed of rigid bodies, constraints,
    %   and joints.
    
    properties
        myBodies; % Structure containing all of the bodies in the system
        myConstraints; % Structure containing all of the constraints in the system
        myNumKinematicConstraints; % Number of kinematic constraints in the system
        myNumDrivingConstraints; % Number of driving constraints in the system
        myBodyIsGround; % Flag indicating if one of the bodies in the system is the ground.
        myR; % Position of all bodies in system
        myRDot; % Velocity of all bodies in system
        myRDDot; % Acceleration of all bodies in system
        myP; % Euler parameters of all bodies in system
        myPDot; % First time derivative of euler parameters of all bodies in system
        myPDDot; % Second time derivative of euler parameters of all bodies in system
        myJMatrixTotal; % Total J matrix. Used for dynamics and inverse dynamics analysis.
        myPhiK; % Vector of kinematic constraints
        myPhiD; % Vector of driving constraints
        myPhiP; % Vector of Euler parameter normalization constraints
        myPhiFull; % Full constraint matrix. Made up of combination of kinematic, driving, and Euler parameter normalization constraints.
        myPhiFullJacobian; % Jacobian of full constraint matrix.
        myPhiPartialR; % Partial derivative of constraint matrix w.r.t. position (not including Euler normalization constraints) 
        myPhiPartialP; % Partial derivative of constraint matrix w.r.t euler parameters (not including Euler normalization constraints)
        myNu; % RHS of velocity equation for entire multibody system
        myGamma; % RHS of acceleration equation for entire multibody system
        myTime; % Current time within system.
        myLagrangeMultipliers; % Lagrange multipliers for constraints.
        myPMatrix; % Euler parametrization constraint matrix
        myInvDynRHS; % RHS of inverse dynamics equations
        myInvDynMatrix; % Matrix that contains partial derivatives of constraints. This is the A matrix when solving Ax = b in inverse dynamics analysis.
    end
    
    properties (Dependent)
        myNumBodies; % Number of bodies in the system
        myNumConstraints; % Number of kinematic and driving constraints in the system
        myMassMatrix; % Mass matrix for the entire system.
        myJMatrix; % Polar moment of inertia matrix for the entire system.
    end
    
    methods
        function obj = multibodySystem()
        end
        function obj = addBody(obj,bodyNumber, bodyType, isGround, mass, bodyLength, JMatrix)
            % mass : double
            %   Mass of body in kg
            %
            % bodyLength : double
            %   Length of body in meters
            
            % Set system level flag indicating that one of the bodies in
            % the system is the ground.
            % NOTE: For ease of populating the Jacobian, make it convention
            % that, if the ground is a body, it is body 1.
            if (isGround == 1)
                if (bodyNumber == 1)
                    obj.myBodyIsGround = 1;
                else
                    disp('ERROR: The ground must be body 1. Please correct this.');
                    clear isGround;
                end
            end
            
            % Add a body to the system
            newBody = body(bodyNumber, bodyType, isGround, mass, bodyLength, JMatrix);
            obj.myBodies{bodyNumber} = newBody;
        end
        function obj = addForce(obj, bodyNumber, force, sBar, forceName)
            % Add a force to a body in the multibody system
            % 
            % Function inputs:
            % bodyNumber : int
            %   Number of the body that you want to apply the force to.
            %            
            % force : 3x1 double
            %   x,y, and z components of force.
            %
            % sBar : 3x1 double
            %   Point of application of the force in the body reference
            %   frame. If the force is applied to the center of mass 
            %   sBar = [0 0 0]'.
            %
            % forceName : string
            %   Name of this force. Optional input.
            %   
            sBar = sBar(:);
            force = force(:);
            nForces = obj.myBodies{bodyNumber}.myNumForces;

            if vargin < 5
                forceName = ['Force ' num2str(nForces+1) ' on body ' num2str(bodyNumber) ];
            end
            
            obj.myBodies{bodyNumber}.addForce(force, sBar, forceName);
        end
        function obj = addTorque(obj, bodyNumber, torque, torqueName)
            % Add a torque to a specific body in the multibody system
            %
            % Function inputs:
            % bodyNumber : int
            %   Body to which you want to apply this torque.
            %
            % torque : 3x1 double
            %   x, y, and z-components of torque.
            %
            % torqueName : string
            %   Name of this torque. Optional input.
            torque = torque(:);
            nTorques = obj.myBodies{bodyNumber}.myNumTorques;

            if vargin < 4
                torqueName = ['Torque ' num2str(nTorques+1) ' on body ' num2str(bodyNumber) ];
            end
            
            obj.myBodies{bodyNumber}.addTorque(torque, torqueName);          
        end
        function obj = kinematicsAnalysis(obj, startTime, endTime, timeStep, displayFlag)
            % Perform kinematics analysis with this multibody system for
            % the given time interval.
            %
            % Function inputs:
            % startTime : double
            %   Start time for kinematics analysis in seconds.
            %
            % endTime : double
            %   End time for kinematics analysis in seconds.
            %
            % timeStep : double
            %   Desired time step for kinematics analysis in seconds.
            
            time = startTime:timeStep:endTime;
            
            for iT = 1:length(time)
                t = time(iT);
                
                % If it is not the first time-step, we must use
                % Newton-Raphson to solve for q. If it is the first
                % time-step, q is given as an initial condition.
                if (t ~= startTime)
                    obj.computeQ(t);
                end
                
                % Check Jacobian at this point. Throw a warning if it is
                % close to singular
                %                 phiFullJacobian = obj.myPhiFullJacobian;
                %                 if norm(phiFullJacobian) < 10^-12
                %                     disp('WARNING: Solution approaching singular configuration!!!');
                %                 end
                
                % Compute qDot
                obj.computeQDot(t);
                
                % Compute gamma and qDDot
                obj.computeQDDot(t);
                
                % Store the position and orientation information for the
                % current time step
                obj.storeSystemState();
                
                if (displayFlag == 1)
                    disp(['Kinematics analysis completed for t = ' num2str(t) ' sec.']);
                end
            end
        end
        function obj = inverseDynamicsAnalysis(obj, startTime, endTime, timestep)
            % Perform an inverse dynamics analysis.
            %
            % Function inputs:
            % startTime : double
            %   Start time for analysis.
            %
            % endTime : double
            %   End time for analysis
            %
            % timestep : double
            %   Time step for analysis.
            time = startTime:timestep:endTime;
            
            for iT = 1:length(time)
                t = time(iT);
                
                % Perform kinematics analysis to get acceleration
                % If it is not the first time-step, we must use
                % Newton-Raphson to solve for q. If it is the first
                % time-step, q is given as an initial condition.
                if (t ~= startTime)
                    obj.computeQ(t);
                end
                
                % Compute qDot
                obj.computeQDot(t);
                
                % Compute gamma and qDDot
                obj.computeQDDot(t);
                
                % Store the position and orientation information for the
                % current time step
                obj.storeSystemState();
                
                % Compute the Lagrange multipliers using the equations of
                % motion
                obj.computeLagrangeMultipliers();
                
            end
               
            
        end
        function obj = computeLagrangeMultipliers(obj)
            % Compute the Lagrange multipliers that will be used for
            % computing the reaction force and torques
            %
            
            % Compute the RHS for linear system of equations for inverse dynamics analysis
            
            
            % Compute the constraint partial derivative matrix
            
            % Compute Lagrange multipliers
            
        end
        function  obj = computeJMatrixTotal(obj)
            % Compute the total J matrix for all bodies.
            
            % Determine number of bodies in system and see the J matrix
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            JMatrixTotal = zeros(4*nBodies,4*nBodies);
            
            % Loop through all bodies. Extract J for the body. Put J into
            % the total matrix.
            iC = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeJpMatrix();
                    JMatrix = obj
                    iC = iC + 1;
                end
                
                
            end
            
            
        end
        function obj = storeSystemState(obj)
            % Store the current state of each body
            
            % Extract matrices that contain the state info for all bodies.
            r = obj.myR;
            rDot = obj.myRDot;
            rDDot = obj.myRDDot;
            p = obj.myP;
            pDot = obj.myPDot;
            pDDot = obj.myPDDot;
            time = obj.myTime;
            
            % Loop through each body. Store the info for each body within
            % its data structure. The resulting data structures will have
            % time and state info across all time steps with each column
            % representing a different time step
            nBodies = obj.myNumBodies;
            for iB = 1:nBodies
                nTimeSteps = obj.myBodies{iB}.myNumTimeSteps;
                obj.myBodies{iB}.myTimeTotal(nTimeSteps + 1) = time;
                obj.myBodies{iB}.myRTotal(:,nTimeSteps + 1) = r(:,iB);
                obj.myBodies{iB}.myRDotTotal(:,nTimeSteps + 1) = rDot(:,iB);
                obj.myBodies{iB}.myRDDotTotal(:,nTimeSteps + 1) = rDDot(:,iB);
                obj.myBodies{iB}.myPTotal(:,nTimeSteps + 1) = p(:,iB);
                obj.myBodies{iB}.myPDotTotal(:,nTimeSteps + 1) = pDot(:,iB);
                obj.myBodies{iB}.myPDDotTotal(:,nTimeSteps + 1) = pDDot(:,iB);
            end
        end
        function [phiFull, phiFullJacobian] = NRfunction(obj, q)
            % Update the state of the system based off the q that is sent
            % in
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                rTemp = q(1:3*nBodies);
                pTemp = q((3*nBodies + 1):7*nBodies);
                
                rGround = [0 0 0]';
                pGround = [1 0 0 0]';
                
                rNew = [rGround; rTemp];
                pNew = [pGround; pTemp];
            else
                nBodies = obj.myNumBodies;
                rNew = q(1:3*nBodies);
                pNew = q((3*nBodies + 1):7*nBodies);
            end
            
            % Reshape into a matrix
            rMatrix = reshape(rNew,[3 obj.myNumBodies]);
            pMatrix = reshape(pNew, [4 obj.myNumBodies]);
            
            % Update the position.
            obj.updateSystemState(rMatrix, [], [], pMatrix, [], [], obj.myTime);
            
            % Now that the system state is updated, recompute phi and the
            % jacobian of phi
            obj.computePhiFull();
            obj.computePhiFullJacobian();
            phiFull = obj.myPhiFull;
            phiFullJacobian = obj.myPhiFullJacobian;
        end
        function obj = computeQ(obj,time)
            % Compute q for all bodies in the system using the
            % Newton-Raphson method. This function must be used at all time
            % steps after the first time step.
            
            % Start with intial guess of q
            rInit = obj.myR;
            pInit = obj.myP;
            
            % Reshape the intial guess of q so it is a vector and remove
            % the body that represents the ground, if one exists.
            if (obj.myBodyIsGround == 1)
                % By convention, body 1 is always the ground, if there is a
                % ground.
                rInit(:,1) = [];
                pInit(:,1) = [];
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            rInitVec = reshape(rInit,[3*nBodies 1]);
            pInitVec = reshape(pInit,[4*nBodies 1]);
            qGuess = [rInitVec; pInitVec];
            
            
            % Use Newton-Raphson method to compute q
            maxIter = 50;
            tol = 1e-9;
            
            iter = 1;
            while iter < maxIter
                % Compute phiFullJacobian and phiFull
                obj.computePhiFull();
                obj.computePhiFullJacobian();
                
                phiFull = obj.myPhiFull;
                phiFullJacobian = obj.myPhiFullJacobian;
                
                % Compute the correction factor
                correction = phiFullJacobian\phiFull;
                
                % Update guess based off correction factor
                qNew = qGuess - correction;
                qGuess = qNew;
                
                % Update the system state so that the new position is qNew.
                % If one of the bodies is the
                % ground you need to put r and p for the ground back into
                % the results, just so we maintain the correct number of total
                % bodies in the system.
                if (obj.myBodyIsGround == 1)
                    nBodies = obj.myNumBodies - 1;
                    rTemp = qNew(1:3*nBodies);
                    pTemp = qNew((3*nBodies + 1):7*nBodies);
                    
                    rGround = [0 0 0]';
                    pGround = [1 0 0 0]';
                    
                    rNew = [rGround; rTemp];
                    pNew = [pGround; pTemp];
                else
                    nBodies = obj.myNumBodies;
                    rNew = qNew(1:3*nBodies);
                    pNew = qNew((3*nBodies + 1):7*nBodies);
                end
                
                % Reshape into a matrix
                rMatrix = reshape(rNew,[3 obj.myNumBodies]);
                pMatrix = reshape(pNew, [4 obj.myNumBodies]);
                
                % Update the position. If this is the last time through the
                % loop, this will allow us to update the system to the
                % final solution.
                obj.updateSystemState(rMatrix, [], [], pMatrix, [], [], time);
                
                % Check to see if the tolerance for convergence has been
                % reached. If it has, break from this loop.
                if (norm(correction) < tol)
                    break
                end
                
                % Update iteration count
                iter = iter + 1;
            end
            
            % Check to see if theta is satisfied.......
            
            
            %             % Perform newton-raphson using fsolve
            %             opt = optimset('fsolve');
            %             opt.GradObj = 'on';
            %             qFinal = fsolve(@(q)obj.NRfunction(q),qGuess);
            %
            %             % Update the system state with the final position
            %             % If one of the bodies is the
            %             % ground you need to put r and p for the ground back into
            %             % the results, just so we maintain the correct number of total
            %             % bodies in the system.
            %             if (obj.myBodyIsGround == 1)
            %                 nBodies = obj.myNumBodies - 1;
            %                 rTemp = qFinal(1:3*nBodies);
            %                 pTemp = qFinal((3*nBodies + 1):7*nBodies);
            %
            %                 rGround = [0 0 0]';
            %                 pGround = [1 0 0 0]';
            %
            %                 rNew = [rGround; rTemp];
            %                 pNew = [pGround; pTemp];
            %             else
            %                 nBodies = obj.myNumBodies;
            %                 rNew = qFinal(1:3*nBodies);
            %                 pNew = qFinal((3*nBodies + 1):7*nBodies);
            %             end
            %
            %             % Reshape into a matrix
            %             rMatrix = reshape(rNew,[3 obj.myNumBodies]);
            %             pMatrix = reshape(pNew, [4 obj.myNumBodies]);
            %
            %             % Update the position. If this is the last time through the
            %             % loop, this will allow us to update the system to the
            %             % final solution.
            %             obj.updateSystemState(rMatrix, [], [], pMatrix, [], [], time);
            %
            %                 %%%%%
            
        end
        
        function obj = computeQDot(obj, time)
            % Compute the velocity and time derivative of the Euler
            % parameters for all bodies in the system.
            obj.myTime = time;
            
            % Compute RHS of velocity equation
            obj.computeNu();
            nu = obj.myNu;
            
            % Compute Jacobian of phi
            obj.computePhiFullJacobian();
            phiFullJacobian = obj.myPhiFullJacobian();
            
            % Solve for qDot
            qDot = phiFullJacobian\nu;
            
            % Parse qDot into rDot and pDot. If one of the bodies is the
            % ground you need to put rDot and pDot for the ground back into
            % the results, just so we maintain the correct number of total
            % bodies in the system.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                rDotTemp = qDot(1:3*nBodies);
                pDotTemp = qDot((3*nBodies + 1):7*nBodies);
                
                rDotGround = [0 0 0]';
                pDotGround = [0 0 0 0]';
                
                rDot = [rDotGround; rDotTemp];
                pDot = [pDotGround; pDotTemp];
            else
                nBodies = obj.myNumBodies;
                rDot = qDot(1:3*nBodies);
                pDot = qDot((3*nBodies + 1):7*nBodies);
            end
            
            % Reshape into a matrix
            rDotMatrix = reshape(rDot,[3 obj.myNumBodies]);
            pDotMatrix = reshape(pDot, [4 obj.myNumBodies]);
            
            % Update the velocities
            obj.updateSystemState([], rDotMatrix, [], [], pDotMatrix, [], time);
        end
        function obj = computeQDDot(obj, time)
            % Compute acceleration and second time derivative of Euler
            % parameters for all bodies in system
            
            % Store time
            obj.myTime = time;
            
            % Compute RHS of acceleration equation
            obj.computeGamma();
            gamma = obj.myGamma;
            
            % Do not waste time computing Jacobian of phi again. It was
            % already computed for the velocity analysis so just grab the
            % previously computed Jacobian.
            phiFullJacobian = obj.myPhiFullJacobian();
            
            % Solve for qDDot
            qDDot = phiFullJacobian\gamma;
            
            % Parse qDDot into rDDot and pDDot. If one of the bodies is the
            % ground you need to put rDDot and pDDot for the ground back into
            % the results, just so we maintain the correct number of total
            % bodies in the system.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                rDDotTemp = qDDot(1:3*nBodies);
                pDDotTemp = qDDot((3*nBodies + 1):7*nBodies);
                
                rDDotGround = [0 0 0]';
                pDDotGround = [0 0 0 0]';
                
                rDDot = [rDDotGround; rDDotTemp];
                pDDot = [pDDotGround; pDDotTemp];
            else
                nBodies = obj.myNumBodies;
                rDDot = qDDot(1:3*nBodies);
                pDDot = qDDot((3*nBodies + 1):7*nBodies);
            end
            
            
            % Reshape into a matrix
            rDDotMatrix = reshape(rDDot,[3 obj.myNumBodies]);
            pDDotMatrix = reshape(pDDot, [4 obj.myNumBodies]);
            
            % Update the velocities
            obj.updateSystemState([], [], rDDotMatrix, [], [], pDDotMatrix, time);
        end
        function obj = computePhiFull(obj)
            % Compute each component of phiFull
            obj.computePhiK();
            obj.computePhiD();
            obj.computePhiP();
            
            phiFull = [obj.myPhiK;
                obj.myPhiD;
                obj.myPhiP];
            obj.myPhiFull = phiFull;
        end
        function obj = computePhiK(obj)
            % Loop through each constraint and check if it is a kinematic
            % constraint. If it is, compute phi for the constraint and add
            % it to the constraint matrix
            
            % Seed phiK
            phiK = zeros(obj.myNumKinematicConstraints,1);
            
            nConst = obj.myNumConstraints;
            
            for iC = 1:nConst
                if (obj.myConstraints{iC}.myIsKinematic == 1)
                    time = obj.myTime;
                    phiFlag = 1;
                    obj.computeConstraintProperties(iC, time, phiFlag, 0, 0, 0, 0);
                    phiK(iC,:) = obj.myConstraints{iC}.myPhi;
                end
            end
            
            % Update myPhiK
            obj.myPhiK = phiK;
        end
        function obj = computePhiD(obj)
            % Loop through each constraint and check if it is a driving
            % constraint. If it is, compute phi for the constraint and add
            % it to the constraint matrix
            
            % Seed phiD
            phiD = zeros(obj.myNumDrivingConstraints,1);
            
            nConst = obj.myNumConstraints;
            
            iD = 1;
            for iC = 1:nConst
                if (obj.myConstraints{iC}.myIsKinematic == 0)
                    time = obj.myTime;
                    phiFlag = 1;
                    obj.computeConstraintProperties(iC, time, phiFlag, 0, 0, 0, 0);
                    phiD(iD,:) = obj.myConstraints{iC}.myPhi;
                    iD = iD + 1;
                end
            end
            
            % Update myPhiD
            obj.myPhiD = phiD;
        end
        function obj = computePhiP(obj)
            % There will be one Euler parameter normalization constraint
            % for each body in the system.
            nBodies = obj.myNumBodies;
            
            % Seed phiP. If one of the bodies in the system is the ground,
            % then phiP has one less row.
            if (obj.myBodyIsGround == 1)
                nPconst = nBodies - 1;
                phiP = zeros(nPconst,1);
            else
                nPconst = nBodies;
                phiP = zeros(nPconst,1);
            end
            
            % Loop through each body and add the value of its Euler
            % parameter normalization constrain (which should be 0) to the
            % overall constraint vector. If the current body is the ground,
            % do not include it.
            iP = 1;
            for iB = 1:nBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    p = obj.myBodies{iB}.myP;
                    phiP(iP,:) = p'*p - 1;
                    iP = iP + 1;
                end
            end
            
            % Update myPhiP
            obj.myPhiP = phiP;
        end
        function obj = computePhiFullJacobian(obj)
            % Compute the jacobian of phiFull (i.e. partial derivative of
            % phiFull w.r.t. the generalized coordinates)
            
            % Seed the Jacobian. There will be 7 less columns and one less in the
            % Jacobian if one of the bodies is the ground.
            nBodies = obj.myNumBodies;
            nKDconst = obj.myNumConstraints;
            
            if (obj.myBodyIsGround == 1)
                phiFullPartialR = zeros((nKDconst + nBodies - 1), 3*(nBodies - 1));
                phiFullPartialP = zeros((nKDconst + nBodies - 1), 4*(nBodies - 1));
            else
                phiFullPartialR = zeros((nKDconst + nBodies), 3*nBodies);
                phiFullPartialP = zeros((nKDconst + nBodies), 4*nBodies);
            end
            
            % Loop through each kinematic and driving constraint. Compute
            % phiPartialR and phiPartialP for each. Insert these the
            % correct location in the Jacobian. This location depends on
            % which bodies are involved in the constraint. If one of the
            % bodies is the ground, then it will not appear in the
            % Jacobian.
            for iC = 1:nKDconst
                % Extract the bodies for this constraint
                bodyI = obj.myConstraints{iC}.myBodyI;
                bodyJ = obj.myConstraints{iC}.myBodyJ;
                
                % Compute phiPartialR and phiPartialP
                time = obj.myTime;
                phiPartialRFlag = 1;
                phiPartialPFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, 0, phiPartialRFlag, phiPartialPFlag);
                phiPartialR = obj.myConstraints{iC}.myPhiPartialR;
                phiPartialP = obj.myConstraints{iC}.myPhiPartialP;
                
                % If one of the bodies in the system is the ground adjust
                % the body numbers accordingly.
                if (obj.myBodyIsGround == 1)
                    if (obj.myBodies{bodyI}.myIsGround == 1)
                        % Only care about bodyJ in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. This only works because I
                        % am forcing the ground to be body 1.
                        phiFullPartialR(iC,(3*(bodyJ-1) - 2):3*(bodyJ-1)) = phiPartialR;
                        phiFullPartialP(iC,(4*(bodyJ-1) - 3):4*(bodyJ-1)) = phiPartialP;
                        
                    elseif (obj.myBodies{bodyJ}.myIsGround == 1)
                        % Only care about bodyI in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian.
                        phiFullPartialR(iC,(3*(bodyI-1) - 2):3*(bodyI-1)) = phiPartialR;
                        phiFullPartialP(iC,(4*(bodyI-1) - 3):4*(bodyI-1)) = phiPartialP;
                        
                    else
                        % Include both bodyI and bodyJ in this case because
                        % neither is the ground, but still decrease body
                        % number by 1.
                        % Populate the partial derivative w.r.t. position.
                        phiFullPartialR(iC,(3*(bodyI-1) - 2):3*(bodyI-1)) = phiPartialR(1:3);
                        phiFullPartialR(iC,(3*(bodyJ-1) - 2):3*(bodyJ-1)) = phiPartialR(4:6);
                        
                        % Populate the partial derivative w.r.t. orientation.
                        phiFullPartialP(iC,(4*(bodyI-1) - 3):4*(bodyI-1)) = phiPartialP(1:4);
                        phiFullPartialP(iC,(4*(bodyJ-1) - 3):4*(bodyJ-1)) = phiPartialP(5:8);
                        
                    end
                else
                    % Populate the partial derivative w.r.t. position.
                    phiFullPartialR(iC,(3*bodyI - 2):3*bodyI) = phiPartialR(1:3);
                    phiFullPartialR(iC,(3*bodyJ - 2):3*bodyJ) = phiPartialR(4:6);
                    
                    % Populate the partial derivative w.r.t. orientation.
                    phiFullPartialP(iC,(4*bodyI - 3):4*bodyI) = phiPartialP(1:4);
                    phiFullPartialP(iC,(4*bodyJ - 3):4*bodyJ) = phiPartialP(5:8);
                end
            end
            
            % Loop through each body and finish populating phiFullPartialP
            % by computing the partial derivative for each of the Euler
            % parameter normalization constraints. There is no need to do
            % anything else with phiFullPartialR because the partial
            % derivative of the Euler parameter constraints are all 0 w.r.t
            % R.
            iP = 1;
            for iB = 1:nBodies
                p = obj.myBodies{iB}.myP;
                if (obj.myBodyIsGround == 1)
                    if (obj.myBodies{iB}.myIsGround == 0)
                        phiFullPartialP((nKDconst + iP),(4*iP - 3):(4*iP)) = 2*p';
                        iP = iP + 1;
                    end
                else
                    phiFullPartialP((nKDconst + iP),(4*iP - 3):(4*iP)) = 2*p';
                    iP = iP + 1;
                end
            end
            
            % Append phiFullPartialR and phiFullPartialP together
            phiFullJacobian = [phiFullPartialR phiFullPartialP];
            obj.myPhiFullJacobian = phiFullJacobian;
            
        end
        function obj = computeNu(obj)
            % Compute RHS of velocity equation for the entire system.
            
            % Seed nu. The number of inputs into nu should be equal to the
            % number of total constraints (kinematic + driving + Euler
            % parameter normalization)
            nKDconst = obj.myNumConstraints;
            nBodies = obj.myNumBodies;
            if (obj.myBodyIsGround == 1) % If one of the bodies is the ground, it does not contibute to nu
                nuLength = nKDconst + nBodies - 1;
            else
                nuLength = nKDconst + nBodies;
            end
            nuTotal = zeros(nuLength,1);
            
            % Loop through each kinematic and driving constraint and compute nu.
            % The rest of the inputs into nu (for the Euler parameter
            % normalization constraints) should remain zero.
            time = obj.myTime;
            for iC = 1:nKDconst
                nuFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, nuFlag, 0, 0, 0);
                nuTotal(iC,:) = obj.myConstraints{iC}.myNu;
            end
            
            % Update system level Nu
            obj.myNu = nuTotal;
        end
        function obj = computeGamma(obj)
            % Compute RHS of acceleration equation for entire multibody
            % system
            
            time = obj.myTime;
            % Seed gamma. The length of gamma is equal to the total number
            % of constraints (kinematic + driving + Euler param norm).
            % Remember, the ground does not contribute to the Euler param
            % norm constraint.
            nKDconst = obj.myNumConstraints;
            nBodies = obj.myNumBodies;
            if (obj.myBodyIsGround == 1)
                gammaLength = nKDconst + nBodies - 1;
            else
                gammaLength = nKDconst + nBodies;
            end
            gammaTotal = zeros(gammaLength,1);
            
            % Loop through each constraint and comute gamma
            for iC = 1:nKDconst
                gammaFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, gammaFlag, 0, 0);
                gammaTotal(iC,:) = obj.myConstraints{iC}.myGamma;
            end
            
            % Loop through each body and compute gamma. If the body is the
            % ground, skip it.
            if (obj.myBodyIsGround == 1)
                iP = 1;
                for iB = 1:nBodies
                    if (obj.myBodies{iB}.myIsGround == 0)
                        pDot = obj.myBodies{iB}.myPDot;
                        gammaTotal((nKDconst + iP),:) = -2*(pDot'*pDot);
                        iP = iP + 1;
                    end
                end
            else
                for iB = 1:nBodies
                    pDot = obj.myBodies{iB}.myPDot;
                    gammaTotal((nKDconst + iB),:) = -2*(pDot'*pDot);
                end
            end
            
            % Update system level gamma
            obj.myGamma = gammaTotal;
        end
        function plot(obj,vargin)
            % Override Matlab's plot command to plot the multibodySystem
            plot.plotSystem(obj,vargin)
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
        function obj = updateSystemState(obj, rMatrix, rDotMatrix, rDDotMatrix, pMatrix, pDotMatrix, pDDotMatrix, time)
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
            
            % Update system variables if they are not empty. If they are
            % empty it means those variables are not being updated at this
            % time.
            if ~isempty(rMatrix)
                obj.myR = rMatrix;
            end
            if ~isempty(rDotMatrix)
                obj.myRDot = rDotMatrix;
            end
            if ~isempty(rDDotMatrix)
                obj.myRDDot = rDDotMatrix;
            end
            if ~isempty(pMatrix)
                obj.myP = pMatrix;
            end
            if ~isempty(pDotMatrix)
                obj.myPDot = pDotMatrix;
            end
            if ~isempty(pDDotMatrix)
                obj.myPDDot = pDDotMatrix;
            end
            obj.myTime = time;
            
            % Update individual bodies. If the desired matrix is empty,
            % this means that specific value is not being updated at the
            % moment. Just keep this values the same.
            nBodies = obj.myNumBodies;
            for iB = 1:nBodies
                if ~isempty(pMatrix)
                    p = pMatrix(:,iB);
                else
                    p = obj.myBodies{iB}.myP;
                end
                if ~isempty(pDotMatrix)
                    pDot = pDotMatrix(:,iB);
                else
                    pDot = obj.myBodies{iB}.myPDot;
                end
                if ~isempty(pDDotMatrix)
                    pDDot = pDDotMatrix(:,iB);
                else
                    pDDot = obj.myBodies{iB}.myPDDot;
                end
                if ~isempty(rMatrix)
                    r = rMatrix(:,iB);
                else
                    r = obj.myBodies{iB}.myR;
                end
                if ~isempty(rDotMatrix)
                    rDot = rDotMatrix(:,iB);
                else
                    rDot = obj.myBodies{iB}.myRDot;
                end
                if ~isempty(rDDotMatrix)
                    rDDot = rDDotMatrix(:,iB);
                else
                    rDDot = obj.myBodies{iB}.myRDDot;
                end
                obj.myBodies{iB}.updateBody(p, pDot, pDDot, r, rDot, rDDot, time);
            end
        end
        function obj = addBasicConstraint(obj,isKinematic,constraintType,attributes)
            % Add a constraint to the multibody system
            %
            % Function inputs:
            % isKinematic : int
            %   Flag to determine if this is a kinematic or driving
            %   constraint. 1 if constraint is kinematic. 0 if constraint
            %   is driving
            %
            % constraintType : int
            %   Type of basic constraint. Possible options 'dp1', 'dp2',
            %   'd', or 'cd'
            %
            % attributes : struct
            %   Structure containing all of the necessary attributes for
            %   your desired constraint. This function will check to make
            %   sure all necessary attributes are provided.
            switch constraintType
                case 'dp1'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'aBarI'} {'aBarJ'} {'ft'} {'ftDot'} {'ftDDot'}];
                    
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            disp(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
                    
                    % Add a DP1 constraint between two bodies in the system
                    a = attributes;
                    newConstraint = DP1constraint(a.bodyI, a.bodyJ, a.aBarI, a.aBarJ, a.ft, a.ftDot, a.ftDDot, a.constraintName);
                    
                case 'dp2'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'aBarI'} {'sBarIP'} {'sBarJQ'} {'ft'} {'ftDot'} {'ftDDot'}];
                    
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            disp(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
                    
                    % Add a DP2 constraint between two bodies in the system
                    a = attributes;
                    newConstraint = DP2constraint(a.bodyI, a.bodyJ, a.aBarI, a.sBarIP, a.sBarJQ, a.ft, a.ftDot, a.ftDDot, a.constraintName);
                    
                case 'd'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'sBarIP'} {'sBarJQ'} {'ft'} {'ftDot'} {'ftDDot'}];
                    
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            disp(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
                    
                    % Add a D constraint between two bodies in the system
                    a = attributes;
                    newConstraint = Dconstraint(a.bodyI, a.bodyJ, a.sBarIP, a.sBarJQ, a.ft, a.ftDot, a.ftDDot, a.constraintName);
                    
                case 'cd'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'coordVec'} {'sBarIP'} {'sBarJQ'} {'ft'} {'ftDot'} {'ftDDot'}];
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            disp(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
                    
                    % Add a CD constraint between two bodies in the system
                    a = attributes;
                    newConstraint = CDconstraint(a.bodyI, a.bodyJ, a.coordVec, a.sBarIP, a.sBarJQ, a.ft, a.ftDot, a.ftDDot, a.constraintName);
                    
            end
            
            % Set flag for kinematic vs driving constraint
            newConstraint.myIsKinematic = isKinematic;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if (isKinematic)
                if isempty(obj.myNumKinematicConstraints)
                    obj.myNumKinematicConstraints = 0;
                end
                obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 1;
            else
                if isempty(obj.myNumDrivingConstraints)
                    obj.myNumDrivingConstraints = 0;
                end
                obj.myNumDrivingConstraints = obj.myNumDrivingConstraints + 1;
            end
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with new constraint
            obj.myConstraints{nConst + 1} = newConstraint;
            
            
        end
        function obj = computeConstraintProperties(obj, constraintNumber, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
            % Determine type of constraint
            cType = obj.myConstraints{constraintNumber}.myConstraintType;
            
            % Call the correct constraint class based off constraint type.
            switch cType
                case 'DP1'
                    obj.myConstraints{constraintNumber}.computeDP1constraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
                    
                case 'DP2'
                    obj.myConstraints{constraintNumber}.computeDP2constraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
                    
                case 'D'
                    obj.myConstraints{constraintNumber}.computeDconstraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
                    
                case 'CD'
                    obj.myConstraints{constraintNumber}.computeCDconstraint(obj, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag);
            end
        end
        
    end
    % Methods block with no attributes
    methods
        function myNumBodies = get.myNumBodies(obj)
            % Calculate number of bodies in system
            myNumBodies = length(obj.myBodies);
        end
        function myNumConstraints = get.myNumConstraints(obj)
            % Calculate number of constraints in system
            myNumConstraints = length(obj.myConstraints);
        end
        function myMassMatrix = get.myMassMatrix(obj)
            % Each time a new body is added to the system, update thte mass
            % matrix.
            
            % Seed the mass matrix. Loop through each body and place the
            % mass matrix of each body in the total mass matrix, unless the
            % body is the ground. Then it contributes nothing to the mass
            % matrix.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                if (nBodies == 0)
                    myMassMatrix = [];
                else
                    myMassMatrix = zeros(3*nBodies,3*nBodies);
                    iEntries = 1;
                    for iB = 1:(nBodies+1)
                        if (obj.myBodies{iB}.myIsGround == 0)
                            bodyMass = obj.myBodies{iB}.myMassMatrix;
                            myMassMatrix((3*iEntries - 2):3*iEntries,(3*iEntries - 2):3*iEntries) = bodyMass;
                            iEntries = iEntries + 1;
                        end
                    end
                end
            else
                nBodies = obj.myNumBodies;
                myMassMatrix = zeros(3*nBodies,3*nBodies);
                for iB = 1:nBodies
                    bodyMass = obj.myBodies{iB}.myMassMatrix;
                    myMassMatrix((3*iB-2):3*iB,(3*iB-2):3*iB) = bodyMass;
                end
            end
            
            
        end
        function myJMatrix = get.myJMatrix(obj)
            % Each time a new body is added to the system, update the polar
            % moment of inertia matrx.
            
                        % Seed the J matrix. Loop through each body and place the
            % J matrix of each body in the total J matrix, unless the
            % body is the ground. Then it contributes nothing to the J
            % matrix.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                if (nBodies == 0)
                    myJMatrix = [];
                else
                    myJMatrix = zeros(3*nBodies,3*nBodies);
                    iEntries = 1;
                    for iB = 1:(nBodies+1)
                        if (obj.myBodies{iB}.myIsGround == 0)
                            bodyJMatrix = obj.myBodies{iB}.myJMatrix;
                            myJMatrix((3*iEntries - 2):3*iEntries,(3*iEntries - 2):3*iEntries) = bodyJMatrix;
                            iEntries = iEntries + 1;
                        end
                    end
                end
            else
                nBodies = obj.myNumBodies;
                myJMatrix = zeros(3*nBodies,3*nBodies);
                for iB = 1:nBodies
                    bodyJMatrix = obj.myBodies{iB}.myMassMatrix;
                    myJMatrix((3*iB-2):3*iB,(3*iB-2):3*iB) = bodyJMatrix;
                end
            end
        end
    end
    
end