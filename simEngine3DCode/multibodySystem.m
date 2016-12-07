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
        myJpMatrixTotal; % Total Jp matrix. Used for dynamics and inverse dynamics analysis.
        myForceVector; % Vector of all the forces acting on all the bodies in the system.
        myTorqueHatVector; % Vector of all the torqueHats for all the bodies in the system.
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
        myTimeTotal; % Vector containing total time steps throughout simulation
        myConstraintLagrangeMultipliers; % Lagrange multipliers for constraints.
        myEulerParamLagrangeMultipliers; % Lagrange multipliers for Euler parameter normalization constraints.
        myPMatrix; % Euler parametrization constraint matrix
        myInvDynRHS; % RHS of inverse dynamics equations
        myInvDynMatrix; % Matrix that contains partial derivatives of constraints. This is the A matrix when solving Ax = b in inverse dynamics analysis.
        myInitCondFlag; % Flag indicating if the initial conditions satisfy the level zero and level one constraints.
        myLHSforEOM; % Matrix for the left hand side of the matrix-form of the Newton-Euler equations of motion.
        myRHSforEOM; % Matrix for the right hand side of the matrix-form of the Newton-Euler equations of motion.
        myRTotal; % 3n x M matrix containing the positions throughout a dynamics analysis
                  % n = Number of bodies (including ground)
                  % M = Number of time steps
        myRDotTotal; % 3n x M matrix containing the velocities throughout a dynamics analysis
        myRDDotTotal; % 3n x M matrix containing the accelerations throughout a dynamics analysis
        myPTotal; % 4n x M matrix containing the Euler parameters throughout a dynamics analysis
        myPDotTotal; % 4n x M matrix containing the 1st time derivative of the Euler parameters throughout a dynamics analysis
        myPDDotTotal; % 4n x M matrix containing the 2nd time derivative of the Euler parameters throughout a dynamics analysis
        myGMatrix; % Matrix that is computed at each iteration of the numerical integration. Used to compute the correction factor.
        myPsi; % Iteration matrix for Quasi-Newton method
        myFiniteDiffPsi; % Psi computed using finite differences.
        myIterCount; % Iteration count for current time step.
        myIterCountTotal; % Total number of iterations for all timesteps
        myVelocityConstraintViolationTotal; % nConstraint x nTimeSteps matrix containing the velocity constraint violation for each step of the simulation.
        myConstraintForcePartialR; % Partial derivative of constraint forces w.r.t r
        myConstraintForcePartialP; % Partial derivative of constraint forces w.r.t p
        myConstraintTorquePartialR; % Partial derivative of constraint torques w.r.t r
        myConstraintTorquePartialP; % Partial derivative of constraint torques w.r.t p
        myTorqueHatPartialP; % Partial derivative of torqueHat w.r.t p
        myTorqueHatPartialPDot; % Partial derivative of torqueHat w.r.t pDot       
        myJpDDotPartialP; % Partial derivative of [Jp*pDDot] w.r.t p for all bodies.
        myEulerConstraintPartialP; % Partial derivative of [P'*lambdaP] w.r.t. p.
    end
    
    properties (Dependent)
        myNumBodies; % Number of bodies in the system
        myNumBodiesMinusGround; % Number of bodies in the system, not
        % including the ground. If there is not
        % ground in the system, myNumBodiesMinusGround = myNumBodies.
        myNumConstraints; % Number of kinematic and driving constraints in the system
        myNumTimeSteps; % Number of time steps that have been completed in simulation.
        myMassMatrix; % Mass matrix for the entire system.
        myJMatrix; % Polar moment of inertia matrix for the entire system.
    end
    
    methods
        function obj = multibodySystem()
        end
        function obj = addBody(obj,bodyNumber, bodyType, isGround, mass, bodyLength, JMatrix, gravityDirection)
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
                    error('ERROR: The ground must be body 1. Please correct this.');
                end
            end
            
            % Add a body to the system
            newBody = body(bodyNumber, bodyType, isGround, mass, bodyLength, JMatrix, gravityDirection);
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
            
            if nargin < 5
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
            
            if nargin < 4
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
            
            % Check to make sure this is a fully constrained system
            nConstTotal = obj.myNumConstraints + obj.myNumBodiesMinusGround;
            if (nConstTotal < 7*obj.myNumBodiesMinusGround)
                error('System not fully constrained. Cannot perform inverse dynamics analysis')
            elseif (nConstTotal > 7*obj.myNumBodiesMinusGround)
                error('System over-constrained. Cannot perform inverse dynamics analysis')
            end
            
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
                obj.storeSystemState(iT);
                
                if (displayFlag == 1)
                    disp(['Kinematics analysis completed for t = ' num2str(t) ' sec.']);
                end
            end
        end
        function obj = setInitialPose(obj, rInitial, pInitial, assemblyAnalysisFlag)
            % Set the initial position and orientation for each body in the
            % system. Then, perform an assembly analysis to optimize the
            % inital pose.
            %
            % Function inputs:
            % rInitial: 3*nBodies x 1 vector
            %   Vector of the initial position for each body, including the ground.
            %
            % pInitial: 4*nBodies x 1 vector
            %   Vector of the initial orientation for each body, including
            %   the ground.
            %
            % assemblyAnalysisFlag : int
            %   Flag to set whether or not the user wants to perform an
            %   assembly analysis before simulating.
            t = 0;
            obj.updateSystemState( rInitial, [], [], pInitial, [], [], t);
            
            if (assemblyAnalysisFlag == 1)
                obj.assemblyAnalysis();
            end
        end
        function obj = computeAndSetInitialVelocities(obj, knownBodies, knownInitialRDot, knownInitialPDot)
            % Compute initial velocities for unknown bodies.
            %
            % Function inputs:
            % knownBodies : vector
            %   Vector containing the list of bodies that have a known
            %   initial velocity.
            %
            % knownInitalRDot : 3*length(knownBodies) x 1 vector
            %   Vector containing the initial translation velocities for
            %   the bodies within known initial velocities.
            %
            % knownInitialPDot : 4*length(knownBodies) x 1 vector
            %   Vector containing the initial rotational velocities for
            %   the bodies within known initial velocities.
            
            nBodies = obj.myNumBodiesMinusGround;
            
            % Compute the necessary jacobians and the Pmatrix
            obj.computePhiPartialR();
            phiR = obj.myPhiPartialR;
            
            obj.computePhiPartialP();
            phiP = obj.myPhiPartialP;
            
            obj.computePMatrix();
            Pmatrix = obj.myPMatrix;
            
            obj.computeNu();
            nu = obj.myNu;
            nuConst = nu(1:obj.myNumConstraints); % Only want the portion of nu related to the constraints.
            
            % Remove portions of these matrices that are related to bodies
            % with known velocities
            bodiesVec = 1:obj.myNumBodiesMinusGround;
            if (obj.myBodyIsGround == 1)
                bodiesToRemove = knownBodies(:) - 1; % Subtract 1 because the jacobians and the P matrix do not include the ground.
            else
                bodiesToRemove = knownBodies(:);
            end
            
            bodiesToKeep = [];
            for i = 1:length(bodiesToRemove)
                if (i == 1)
                    bodiesToKeep = find(bodiesToRemove(i) ~= bodiesVec);
                else
                    temp = find(bodiesToRemove(i) ~= bodiesVec);
                    bodiesToKeep = [bodiesToKeep temp];
                end
            end
            
            if isempty(bodiesToRemove)
                phiRUnknown = phiR;
                phiPUnknown = phiP;
                PmatrixUnknown = Pmatrix;
                
                % Create matrix to solve for the unknown velocities
                % Left hand side
                nBodiesNew = nBodies - length(bodiesToRemove);
                LHSmatrix = [phiRUnknown phiPUnknown;
                    zeros(nBodiesNew,3*nBodiesNew) PmatrixUnknown];
                
                % Right hand side
                firstEntry = nuConst;
                secondEntry = zeros(nBodiesNew,1);
                RHS = [firstEntry;
                    secondEntry];
                
            else
                phiRKnown = phiR(:,(3*bodiesToRemove - 2):3*bodiesToRemove);
                phiRUnknown = phiR;
                phiRUnknown(:,(3*bodiesToRemove - 2):3*bodiesToRemove) = [];
                
                phiPKnown = phiP(:,(4*bodiesToRemove - 3):4*bodiesToRemove);
                phiPUnknown = phiP;
                phiPUnknown(:,(4*bodiesToRemove - 3):4*bodiesToRemove) = [];
                
                PmatrixKnown = Pmatrix(bodiesToRemove,(4*bodiesToRemove - 3):4*bodiesToRemove);
                
                PmatrixUnknown = Pmatrix;
                PmatrixUnknown(bodiesToRemove,:) = [];
                PmatrixUnknown(:,(4*bodiesToRemove - 3):4*bodiesToRemove) = [];
                
                % Create matrix to solve for the unknown velocities
                % Left hand side
                nBodiesNew = nBodies - length(bodiesToRemove);
                LHSmatrix = [phiRUnknown phiPUnknown;
                    zeros(nBodiesNew,3*nBodiesNew) PmatrixUnknown];
                
                % Right hand side
                firstEntry = nuConst - phiRKnown*knownInitialRDot - phiPKnown*knownInitialPDot;
                secondEntry = zeros(nBodiesNew,1);
                RHS = [firstEntry;
                    secondEntry];
            end
            
            
            % Solve for the unknown velocities.
            vel = LHSmatrix\RHS;
            
            % Extract translational and rotational velocities
            velR = vel(1:3*nBodiesNew);
            velP = vel((3*nBodiesNew + 1):7*nBodiesNew);
                        
            % Add in the known velocities
            velRTotal = zeros(3*nBodies,1);
            velPTotal = zeros(4*nBodies,1);
            unknownCount = 1;
            knownCount = 1;
            for iB = 1:nBodies
                if any(bodiesToRemove == iB)
                    velRTotal((3*iB - 2):3*iB) = knownInitialRDot((3*knownCount - 2):3*knownCount);
                    velPTotal((4*iB - 3):4*iB) = knownInitialPDot((4*knownCount - 3):4*knownCount);
                    knownCount = knownCount + 1;
                else
                    velRTotal((3*iB - 2):3*iB) = velR((3*unknownCount - 2):3*unknownCount);
                    velPTotal((4*iB - 3):4*iB) = velP((4*unknownCount - 3):4*unknownCount);
                    unknownCount = unknownCount + 1;
                end
            end
                    
            
            % Add in the ground if it is in the system
            if (obj.myBodyIsGround == 1)
                rDotTemp = velRTotal;
                pDotTemp = velPTotal;
                
                rDotInitial = [0; 0; 0; rDotTemp];
                pDotInitial = [0; 0; 0; 0; pDotTemp];
            else
                rDotInitial = velRTotal;
                pDotInitial = velPTotal;
            end
            
             % Update the system state to store these computed velocities.
             time = 0;
             obj.updateSystemState([], rDotInitial, [], [], pDotInitial, [], time);
        end        
        function obj = removeRedundantConstraints(obj)
            % Remove redundant constraints from an over-constrained system.
            % A successful assembly analysis must have been performed prior
            % to running this function.
            
            % Evalute Jacobian of kinematic constraint equations
            nKinematicConst = obj.myNumKinematicConstraints;
            obj.computePhiPartialR();
            phiPartialR = obj.myPhiPartialR(1:nKinematicConst,:);
            
            obj.computePhiPartialP();
            phiPartialP = obj.myPhiPartialP(1:nKinematicConst,:);
            
            phiKJacobian = [phiPartialR phiPartialP];
            
            
            % Determine the rank of this Jacobian
            rankOfJacobian =  rank(phiKJacobian);
            
            A = rref(phiKJacobian);
            
            [L, U, P, Q] = simEngine3DUtilities.gecp(phiKJacobian);
            
        end
        function obj = dynamicsAnalysis(obj, startTime, endTime, timestep, order, iterationMethod, displayFlag)
            % Perform a dynamics analysis.
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
            %
            % order : int
            %   Order of the BDF method to be used, as specificed by the user
            %
            % iterationMethod : string
            %   Method used to compute iteration matrix.
            %   Options are:
            %   'newtonRaphson'
            %   'modifiedNewton'
            %   'quasiNewton'
            %
            % displayFlag : int
            %   Flag for user to indicate if they want to display when each
            %   time step of analysis has been completed.
            
            % Check to make sure the system is not overconstrained.
            nConstTotal = obj.myNumConstraints + obj.myNumBodiesMinusGround;
            if (nConstTotal > 7*obj.myNumBodiesMinusGround)
                disp('System over-constrained. Attempting to remove redundant constraints.')
                obj.removeRedundantConstraints();
            end
                       
            % Compute time vector for analysis
            time = startTime:timestep:endTime;
            
            % Allocate space for matrices that will be populated throughout
            % dynamics analysis
            nB = obj.myNumBodies;
            nC = obj.myNumConstraints;
            nT = length(time);
            obj.myTimeTotal = zeros(1,nT);
            obj.myRTotal = zeros(3*nB,nT);
            obj.myRDotTotal = zeros(3*nB,nT);
            obj.myRDDotTotal = zeros(3*nB,nT);
            obj.myPTotal = zeros(4*nB,nT);
            obj.myPDotTotal = zeros(4*nB,nT);
            obj.myPDDotTotal = zeros(4*nB,nT);
            obj.myIterCountTotal = zeros(1,nT);
            obj.myVelocityConstraintViolationTotal = zeros(nC, nT);
            
            % Allocate space within each body as well.
            for iB = 1:nB
                obj.myBodies{iB}.myTimeTotal = zeros(1,nT);
                obj.myBodies{iB}.myRTotal = zeros(3,nT);
                obj.myBodies{iB}.myRDotTotal = zeros(3,nT);
                obj.myBodies{iB}.myRDDotTotal = zeros(3,nT);
                obj.myBodies{iB}.myPTotal = zeros(4,nT);
                obj.myBodies{iB}.myPDotTotal = zeros(4,nT);
                obj.myBodies{iB}.myPDDotTotal = zeros(4,nT);
                obj.myBodies{iB}.myConstraintForcesTotal = zeros(3*nC,nT);
                obj.myBodies{iB}.myConstraintTorquesTotal = zeros(4*nC,nT);
                obj.myBodies{iB}.myConstraintTorquesOmegaTotal = zeros(3*nC,nT);
            end
            
            for iT = 1:length(time)
                t = time(iT);
                
                if (iT == 1)
                    % For first iteration, make sure the initial conditions
                    % satisfy the level zero and level one constraint
                    % equations.
                    obj.checkInitialConditions();
                    if (obj.myInitCondFlag == 0)
                        disp('Initial conditions do not satisfy constraint equations. Please correct this.');
                        break;
                    end
                    
                    % Also for the first iteration, compute qDDot and the
                    % Lagrange multipliers for the initial position and
                    % velocity conditions.
                    obj.solveFirstStepDynamicsAnalysis(t);
                    
                else
                    % If this is not the first time step use a BDF method
                    % to perform numerical integration and get the
                    % accelerations and lagrange mutlipliers for the first
                    % time step. Within this function we will also compute
                    % the current position and velocity.
                    obj.performNumericalIntegrationDynamicsAnalysis(t, timestep, order, iT, iterationMethod);
                    
                end
                
                % Store total number of iterations for previous time step
                if (iT == 1)
                    obj.myIterCountTotal(iT) = 0;
                else
                    obj.myIterCountTotal(iT) = obj.myIterCount;
                end
                
                % Store the current position, velocity, and acceleration of
                % the body for this time step
                obj.storeSystemState(iT);
                
                % Compute the reaction forces from the Lagrange multipliers
                obj.computeConstraintForces();
                obj.computeConstraintTorques();
                
                % Store the reaction forces
                obj.storeConstraintForcesAndTorques(t,iT);
                
                % Compute the violation of the velocity constraints. No
                % need to recompute phiPartialR and phiPartialP within this
                % function because we computed them to compute the
                % constraint forces and the constraint torques.
                obj.computeVelocityConstraintViolation(iT);
                
                if (displayFlag == 1)
                    disp(['Dynamics analysis completed for t = ' num2str(t) ' sec.']);
                end
            end
        end
        function obj = computeVelocityConstraintViolation(obj, stepNumber)
            % Compute violation of velocity constraint
            %
            % Function inputs:
            % stepNumber : int
            %   Current step number of the dynamics analysis.
            
            % Check the velocity constraint.
            % Compute necessary parameters. We computed phiPartialR and
            % phiPartialP prior to this function so no need to recompute
            % them.
            obj.computeNu();
            
            % Extract necessary parameters.
            phiPartialR = obj.myPhiPartialR;
            phiPartialP = obj.myPhiPartialP;
            nu = obj.myNu;
            rDot = obj.myRDot;
            pDot = obj.myPDot;
            
            % Remove ground from system if it exists
            if (obj.myBodyIsGround == 1)
                % rDot and pDot include the ground so remove them.
                rDotVec = rDot(4:end,1);
                pDotVec = pDot(5:end,1);
            else
                rDotVec = rDot;
                pDotVec = pDot;
            end
            
            % Extract the portion of nu that relates to the constraints
            nuConst = nu(1:obj.myNumConstraints);
            
            % Check to ensure this constraint is satisfied
            velConstViolation = phiPartialR*rDotVec + phiPartialP*pDotVec - nuConst;
            obj.myVelocityConstraintViolationTotal(:,stepNumber) = velConstViolation;
        end
        function obj = performNumericalIntegrationDynamicsAnalysis(obj, time, stepsize, order, stepNumber, iterationMethod)
            % For this current time step, use a numerical integration
            % technique to compute the position, velocity, accelerations,
            % and Lagrange multipliers.
            %
            % Function inputs :
            % time : double
            %   Time for current timestep
            %
            % stepsize : double
            %   Step size being used for current simulation.
            %
            % order : int
            %   Order of the BDF method to be used, as specificed by the user
            %
            % stepNumber : int
            %   Current step number of the dynamics analysis.
            %
            % iterationMethod : string
            %   Method used to compute iteration matrix.
            %   Options are:
            %   'newtonRaphson'
            %   'modifiedNewton'
            %   'quasiNewton'
            %
            
            % Set variables for convergence
            maxIter = 50;
            tolerance = 1e-3;
            
            % Set the intial guess for the accelerations and lagrange
            % multipliers. These will be the state of the system at the
            % end of the previous time step.
            nBodies = obj.myNumBodiesMinusGround;
            nConst = obj.myNumConstraints;
            
            % If a ground is in the system, remove the ground from the
            % guess of rDDot and pDDot
            rDDotGuess = obj.myRDDotTotal(:,(stepNumber-1));
            pDDotGuess = obj.myPDDotTotal(:,(stepNumber-1));
            pLambdaGuess = obj.myEulerParamLagrangeMultipliers;
            constLambdaGuess = obj.myConstraintLagrangeMultipliers;
            
            % Remove the first body if there is a ground in the system
            if (obj.myBodyIsGround == 1)
                rDDotNew = rDDotGuess(4:end,1);
                rDDotGuess = rDDotNew;
                
                pDDotNew = pDDotGuess(5:end,1);
                pDDotGuess = pDDotNew;
            end
            
            % Create the total guess vector
            zGuess = [rDDotGuess;
                pDDotGuess;
                pLambdaGuess;
                constLambdaGuess];
            
            
            % Begin numerical integration
            iter = 1;
            while iter < maxIter
                
                % Update the Lagrange multipliers with current guess for
                % the Lagrange multipliers
                obj.myConstraintLagrangeMultipliers = constLambdaGuess;
                obj.myEulerParamLagrangeMultipliers = pLambdaGuess;
                
                % Compare the quasi-newton iteration matrix and a numerical
                % approximation of the iteration matrix.
                %                 if (iter == 1) && (stepNumber == 1)
                finiteDiffFlag = 0;
                if (iter == 1) && (finiteDiffFlag == 1)
                    
                    obj.computeQuasiNewtonPsi();
                    
                    % Compute a numerical approximation of the Jacobian
                    if (stepNumber <=2)
                        orderTemp = 1;
                    else
                        orderTemp = 2;
                    end
                    
                    tTemp = time - stepsize;
                    obj.computeFiniteDiffPsi(zGuess, stepsize, orderTemp, tTemp);
                end
                
                % Compute position and velocities using BDF method. Within
                % this function, the position, velocities, and
                % accelerations of the system will be updated. Also,
                % compute G matrix
                if (order == 2) && (stepNumber <= 2)
                    orderTemp = 1;
                    obj.BDFmethodIteration(rDDotGuess, pDDotGuess, orderTemp, stepsize, time, stepNumber);
                    obj.computeGMatrix(stepsize,orderTemp);
                else
                    obj.BDFmethodIteration(rDDotGuess, pDDotGuess, order, stepsize, time, stepNumber);
                    obj.computeGMatrix(stepsize,order);
                end
                
                % Extract residual in the system (g matrix)
                gMatrix = obj.myGMatrix;
                
                % Compute iteration matrix based on the desired method. For
                % modified newton and quasi-newton you only compute the
                % iteration matrix on the first iteration.
                switch iterationMethod
                    case 'newtonRaphson'
                        if (order == 2) && (stepNumber <= 2)
                            orderTemp = 1;
                        else
                            orderTemp = 2;
                        end
                        obj.computeNewtonRaphsonPsi(stepsize, orderTemp);
                    case 'modifiedNewton'
                        if (iter == 1)
                            if (order == 2) && (stepNumber <= 2)
                                orderTemp = 1;
                            else
                                orderTemp = 2;
                            end
                            obj.computeNewtonRaphsonPsi(stepsize, orderTemp);
                        end
                    case 'quasiNewton'
                        if (iter == 1)
                            obj.computeQuasiNewtonPsi();
                        end
                end
                psi = obj.myPsi;
                
                % Solve linear system to get the correction to this guess
                correction = psi\-gMatrix;
                
                % Improve the guess
                zNew = zGuess + correction;
                zGuess = zNew;
                
                % Extract each component of zGuess
                rDDotGuess = zGuess(1:3*nBodies,1);
                pDDotGuess = zGuess((3*nBodies+1):(7*nBodies),1);
                pLambdaGuess = zGuess((7*nBodies+1):(8*nBodies),1);
                constLambdaGuess = zGuess((8*nBodies+1):(8*nBodies+nConst),1);
                
                % Check for convergence
                if (norm(correction) < tolerance)
                    break;
                end
                
                iter = iter + 1;
            end
            
            % Store the total number of iterations
            obj.myIterCount = iter;
            
            % Store the final values of the Lagrange multipliers.
            obj.myConstraintLagrangeMultipliers = constLambdaGuess;
            obj.myEulerParamLagrangeMultipliers = pLambdaGuess;
            
            % Compute the position and velocity again from these final values
            % of the acceleration. The system state is updated within this
            % function.
            if (order == 2) && (stepNumber <= 2)
                orderTemp = 1;
                obj.BDFmethodIteration(rDDotGuess, pDDotGuess, orderTemp, stepsize, time, stepNumber);
            else
                obj.BDFmethodIteration(rDDotGuess, pDDotGuess, order, stepsize, time, stepNumber);
            end
        end
        function obj = computeFiniteDiffPsi(obj, zGuess, stepsize, order, time)
            % Compute finite difference approximation of psi.
            %
            % zGuess : (8*nBodies + nConstraints) x 1 vector
            %   Vector containing the current guess for z. For this finite
            %   difference approximation of psi to be accurate, z must be a
            %   "healthy" set of generalized coordinates and constraints.
            
            delta = 10^-3;
            
            finiteDiffPsi = zeros(length(zGuess),length(zGuess));
            
            % Compute G for this current guess
            nBodies = obj.myNumBodiesMinusGround;
            nConst = obj.myNumConstraints;
            rDDotGuess = zGuess(1:3*nBodies,1);
            pDDotGuess = zGuess((3*nBodies+1):(7*nBodies),1);
            pLambdaGuess = zGuess((7*nBodies+1):(8*nBodies),1);
            constLambdaGuess = zGuess((8*nBodies+1):(8*nBodies+nConst),1);
            
            if (obj.myBodyIsGround == 1)
                obj.updateSystemState([],[],[zeros(3,1); rDDotGuess],[],[],[zeros(4,1); pDDotGuess],time);
            else
                obj.updateSystemState([],[],rDDotGuess,[],[],pDDotGuess,time);
            end
            
            obj.myConstraintLagrangeMultipliers = constLambdaGuess;
            obj.myEulerParamLagrangeMultipliers = pLambdaGuess;
            obj.computeGMatrix(stepsize, order);
            G0 = obj.myGMatrix;
            
            % Loop through each value of zGuess, increment it by a delta,
            % save these new values, compute the new value of G, compute
            % sensitivity of G, place this into the finite differences
            % jacobian, and restore the values of zGuess.
            for iZ = 1:length(zGuess)
                zNew = zGuess;
                zNew(iZ) = zNew(iZ) + delta;
                
                % Store these new values within the system.
                rDDotNew = zNew(1:3*nBodies,1);
                pDDotNew = zNew((3*nBodies+1):(7*nBodies),1);
                pLambdaNew = zNew((7*nBodies+1):(8*nBodies),1);
                constLambdaNew= zNew((8*nBodies+1):(8*nBodies+nConst),1);
                
                if (obj.myBodyIsGround == 1)
                    obj.updateSystemState([],[],[zeros(3,1); rDDotNew],[],[],[zeros(4,1); pDDotNew],time);
                else
                    obj.updateSystemState([],[],rDDotNew,[],[],pDDotNew,time);
                end
                obj.myConstraintLagrangeMultipliers = constLambdaNew;
                obj.myEulerParamLagrangeMultipliers = pLambdaNew;
                
                % Compute the new G matrix
                obj.computeGMatrix(stepsize, order);
                Gnew = obj.myGMatrix;
                
                % Compute the finite differences.
                finiteDiff = (Gnew - G0)/delta;
                
                % Place this finite differences into the jacobian
                finiteDiffPsi(:,iZ) = finiteDiff;
            end
            
            % Store psi
            obj.myFiniteDiffPsi = finiteDiffPsi;
            
            % Restore the system state back to the original values of
            % zGuess.
            if (obj.myBodyIsGround == 1)
                obj.updateSystemState([],[],[zeros(3,1); rDDotGuess],[],[],[zeros(4,1); pDDotGuess],time);
            else
                obj.updateSystemState([],[],rDDotGuess,[],[],pDDotGuess,time);
            end
            
        end
        function obj = computeNewtonRaphsonPsi(obj, stepsize, order)
            % Compute the iteration matrix (psi) for the Newton-Raphson
            % method
            %
            % Function inputs:
            % stepsize : double
            %   Step size being used for current simulation.
            %
            % order : int
            %   Order of the BDF method to be used, as specificed by the user
            
            if (order == 1)
                beta0 = 1;
            elseif (order == 2)
                beta0 = 2/3;
            end
            
            % Extract terms for iteration matrix. All of the terms were
            % previously compute when we computed the Gmatrix, so they do
            % not need to be recomputed.
            massMatrix = obj.myMassMatrix;
            Pmatrix = obj.myPMatrix;
            phiPartialR = obj.myPhiPartialR;
            phiPartialP = obj.myPhiPartialP;
            nBodies = obj.myNumBodiesMinusGround;
            nConst = obj.myNumConstraints;
            
            % Compute the partial derivative of the constraint forces  and torques w.r.t
            % r and p
            obj.computePartialDerivativeConstraintForcesAndTorques();
            constraintForcePartialR = obj.myConstraintForcePartialR;
            constraintForcePartialP =  obj.myConstraintForcePartialP;
            constraintTorquePartialR = obj.myConstraintTorquePartialR;
            constraintTorquePartialP = obj.myConstraintTorquePartialP;
            
            % Set the partial derivative of the externally applied
            % forces w.r.t r, rDot, p, and pDot to zero.
            forcePartialR = zeros(3*nBodies, 3*nBodies);
            forcePartialRDot = zeros(3*nBodies, 3*nBodies);
            forcePartialP = zeros(3*nBodies, 4*nBodies);
            forcePartialPDot = zeros(3*nBodies, 4*nBodies);
            
            % Set the partial derivative of the externally applied
            % torques w.r.t r and rDot to zero.
            tauHatPartialR = zeros(4*nBodies, 3*nBodies);
            tauHatPartialRDot = zeros(4*nBodies, 3*nBodies);
            
            % Compute the partial derivative of the externally applied
            % torques w.r.t. p
            obj.computeTorqueHatPartialP();
            tauHatPartialP = obj.myTorqueHatPartialP;
            
            % Compute the partial derivative of the externally applied
            % torques w.r.t. pDot
            obj.computeTorqueHatPartialPDot();
            tauHatPartialPDot = obj.myTorqueHatPartialPDot;

            % Compute partial derivative of Euler normalization constraint
            % w.r.t p.
            obj.computeEulerConstraintPartialP();
            eulerConstPartialP = obj.myEulerConstraintPartialP;
            
            % Compute partial derivative of Jp*pDDot w.r.t p.
            obj.computeJpDDotPartialP();
            JpDDotPartialP = obj.myJpDDotPartialP;
            
            % Compute JpMatrix
            obj.computeJpMatrixTotal();
            JpMatrix = obj.myJpMatrixTotal;
            
            % Compute components of psi
            psi11 = massMatrix + ...
                stepsize^2*beta0^2*constraintForcePartialR - ...
                stepsize^2*beta0^2*forcePartialR - ...
                stepsize*beta0*forcePartialRDot;
            
            psi12 = stepsize^2*beta0^2*constraintForcePartialP - ...
                stepsize^2*beta0^2*forcePartialP - ...
                stepsize*beta0*forcePartialPDot;
            
            psi21 = stepsize^2*beta0^2*constraintTorquePartialR - ...
                stepsize^2*beta0^2*tauHatPartialR - ...
                stepsize*beta0*tauHatPartialRDot;
            
            psi22 = JpMatrix + ...
                stepsize^2*beta0^2*JpDDotPartialP + ...
                stepsize^2*beta0^2*eulerConstPartialP + ...
                stepsize^2*beta0^2*constraintTorquePartialP - ...
                stepsize^2*beta0^2*tauHatPartialP - ...
                stepsize*beta0*tauHatPartialPDot;
            
            % Matrix Assemble!!!
            psi = [psi11, psi12, zeros(3*nBodies,nBodies), phiPartialR';
                psi21, psi22, Pmatrix', phiPartialP';
                zeros(nBodies,3*nBodies), Pmatrix, zeros(nBodies,nBodies), zeros(nBodies,nConst);
                phiPartialR, phiPartialP, zeros(nConst,nBodies), zeros(nConst,nConst)];
            
            obj.myPsi = psi;         
        end
        function obj = computeEulerConstraintPartialP(obj)
            % Compute partial derivative of [P'*lambdaP] w.r.t. p
            lambdaP = obj.myEulerParamLagrangeMultipliers;
            
            % Seed partial derivative matrix
            nBodies = obj.myNumBodiesMinusGround;
            eulerConstPartialP = zeros(4*nBodies,4*nBodies);
            
            % Loop through each body and add it to the partial derivative
            % matrix
            for iB = 1:nBodies
                eulerConstPartialP((4*iB - 3):4*iB,(4*iB - 3):4*iB) = lambdaP(iB)*eye(4);                
            end
            obj.myEulerConstraintPartialP = eulerConstPartialP;
            
        end
        function obj = computeJpDDotPartialP(obj)
             % Compute the partial derivative Jp*pDDot w.r.t. p.
            nBodies = obj.myNumBodiesMinusGround;
            JpDDotPartialPTotal = zeros(4*nBodies, 4*nBodies);
            
            % Loop through all bodies. Compute torqueHatPartial. Add this
            % torqueHat to the overall vector
            iT = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeJpDDotPartialP();
                    JpDDotPartialP = obj.myBodies{iB}.myJpDDotPartialP;
                    JpDDotPartialPTotal((4*iT-3):(4*iT),(4*iT-3):(4*iT)) = JpDDotPartialP;
                    iT = iT + 1;
                end
            end
            obj.myJpDDotPartialP = JpDDotPartialPTotal;
            
        end
        function obj = computeTorqueHatPartialPDot(obj)
            % Compute the partial derivative of the externally applied
            % torque w.r.t. p
            nBodies = obj.myNumBodiesMinusGround;
            torqueHatPartialPDotTotal = zeros(4*nBodies, 4*nBodies);
            
            % Loop through all bodies. Compute torqueHatPartial. Add this
            % torqueHat to the overall vector
            iT = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeTorqueHatPartialPDot();
                    torqueHatPartialPDot = obj.myBodies{iB}.myTorqueHatPartialPDot;
                    torqueHatPartialPDotTotal((4*iT-3):(4*iT),(4*iT-3):(4*iT)) = torqueHatPartialPDot;
                    iT = iT + 1;
                end
            end
            obj.myTorqueHatPartialPDot = torqueHatPartialPDotTotal;
        end
        function obj = computeTorqueHatPartialP(obj)
            % Compute the partial derivative of the externally applied
            % torque w.r.t. p
            nBodies = obj.myNumBodiesMinusGround;
            torqueHatPartialPTotal = zeros(4*nBodies, 4*nBodies);
            
            % Loop through all bodies. Compute torqueHatPartial. Add this
            % torqueHat to the overall vector
            iT = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeTorqueHatPartialP();
                    torqueHatPartialP = obj.myBodies{iB}.myTorqueHatPartialP;
                    torqueHatPartialPTotal((4*iT-3):(4*iT),(4*iT-3):(4*iT)) = torqueHatPartialP;
                    iT = iT + 1;
                end
            end
            obj.myTorqueHatPartialP = torqueHatPartialPTotal;
        end
        function obj = computePartialDerivativeConstraintForcesAndTorques(obj)
            % Compute the partial derivative of constraint forces and
            % torques w.r.t r and p.
            
            % Seed the resulting matrices
            nBodies = obj.myNumBodiesMinusGround;
            nConst = obj.myNumConstraints;
            constForcesTotalPartialR = zeros(3*nBodies,3*nBodies);
            constForcesTotalPartialP = zeros(3*nBodies,4*nBodies);
            constTorquesTotalPartialR = zeros(4*nBodies,3*nBodies);
            constTorquesTotalPartialP = zeros(4*nBodies,4*nBodies);
            
            % Loop through each constraint. Compute the partial derivatives
            % of the reaction forces and torques for this constraint.
            % Place these in the correct location of the overall constraint
            % matrix.
            for iC = 1:nConst
                % Extract the bodies for this constraint
                bodyI = obj.myConstraints{iC}.myBodyI;
                bodyJ = obj.myConstraints{iC}.myBodyJ;
                
                % Compute partial derivatives of for this constraint
                time = obj.myTime;
                constraintForcePartialRFlag = 1;
                constraintForcePartialPFlag = 1;
                constraintTorquePartialRFlag = 1;
                constraintTorquePartialPFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, 0, 0, 0, constraintForcePartialRFlag, constraintForcePartialPFlag, constraintTorquePartialRFlag, constraintTorquePartialPFlag);
                constForcePartialR = obj.myConstraints{iC}.myConstraintForcePartialR;
                constForcePartialP = obj.myConstraints{iC}.myConstraintForcePartialP;
                constTorquePartialR = obj.myConstraints{iC}.myConstraintTorquePartialR;
                constTorquePartialP = obj.myConstraints{iC}.myConstraintTorquePartialP;
                
                % Place partial derivatives into correct section of overall
                % matrix.
                % If one of the bodies in the system is the ground adjust
                % the body numbers accordingly.
                if (obj.myBodyIsGround == 1)
                    bodyJnew = bodyJ - 1;
                    bodyInew = bodyI - 1;
                    if (obj.myBodies{bodyI}.myIsGround == 1)
                        % Only care about bodyJ in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. This only works because I
                        % am forcing the ground to be body 1.
                        constraintForceRJRJ = constForcePartialR;
                        constraintForceRJPJ = constForcePartialP;
                        constraintTorquePJRJ = constTorquePartialR;
                        constraintTorquePJPJ = constTorquePartialP;
                        
                        constForcesTotalPartialR((3*bodyJnew - 2):3*bodyJnew,(3*bodyJnew - 2):3*bodyJnew) = ...
                            constForcesTotalPartialR((3*bodyJnew - 2):3*bodyJnew,(3*bodyJnew - 2):3*bodyJnew) + constraintForceRJRJ;
                        
                        constForcesTotalPartialP((3*bodyJnew - 2):3*bodyJnew,(4*bodyJnew - 3):4*bodyJ) = ...
                            constForcesTotalPartialP((3*bodyJnew - 2):3*bodyJnew,(4*bodyJnew - 3):4*bodyJnew) + constraintForceRJPJ;
                        
                        constTorquesTotalPartialR((4*bodyJnew - 3):4*bodyJnew,(3*bodyJnew - 2):3*bodyJnew) = ...
                            constTorquesTotalPartialR((4*bodyJnew - 3):4*bodyJnew,(3*bodyJnew - 2):3*bodyJnew) + constraintTorquePJRJ;
                        
                        constTorquesTotalPartialP((4*bodyJnew - 3):4*bodyJnew,(4*bodyJnew - 3):4*bodyJnew) = ...
                            constTorquesTotalPartialP((4*bodyJnew - 3):4*bodyJnew,(4*bodyJnew - 3):4*bodyJnew) + constraintTorquePJPJ;
                        
                    elseif (obj.myBodies{bodyJ}.myIsGround == 1)
                        % Only care about bodyI in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian.
                        constraintForceRIRI = constForcePartialR;
                        constraintForceRIPI = constForcePartialP;
                        constraintTorquePIRI = constTorquePartialR;
                        constraintTorquePIPI = constTorquePartialP;
                        
                        constForcesTotalPartialR((3*bodyInew - 2):3*bodyInew,(3*bodyInew - 2):3*bodyInew) = ...
                            constForcesTotalPartialR((3*bodyInew - 2):3*bodyInew,(3*bodyInew - 2):3*bodyInew) + constraintForceRIRI;
                        constForcesTotalPartialP((3*bodyInew - 2):3*bodyInew,(4*bodyInew - 3):4*bodyInew) = ...
                            constForcesTotalPartialP((3*bodyInew - 2):3*bodyInew,(4*bodyInew - 3):4*bodyInew) + constraintForceRIPI;
                        constTorquesTotalPartialR((4*bodyInew - 3):4*bodyInew,(3*bodyInew - 2):3*bodyInew) = ...
                            constTorquesTotalPartialR((4*bodyInew - 3):4*bodyInew,(3*bodyInew - 2):3*bodyInew) + constraintTorquePIRI;
                        constTorquesTotalPartialP((4*bodyInew - 3):4*bodyInew,(4*bodyInew - 3):4*bodyInew) = ...
                            constTorquesTotalPartialP((4*bodyInew - 3):4*bodyInew,(4*bodyInew - 3):4*bodyInew) + constraintTorquePIPI;
                    else
                        % Include both bodyI and bodyJ in this case because
                        % neither is the ground, but still decrease body
                        % number by 1.
                        constraintForceRIRI = constForcePartialR(1:3,1:3);
                        constraintForceRIRJ = constForcePartialR(1:3,4:6);
                        constraintForceRJRI = constForcePartialR(4:6,1:3);
                        constraintForceRJRJ = constForcePartialR(4:6,4:6);
                        
                        constraintForceRIPI = constForcePartialP(1:3,1:4);
                        constraintForceRIPJ = constForcePartialP(1:3,5:8);
                        constraintForceRJPI = constForcePartialP(4:6,1:4);
                        constraintForceRJPJ = constForcePartialP(4:6,5:8);
                        
                        constraintTorquePIRI = constTorquePartialR(1:4,1:3);
                        constraintTorquePIRJ = constTorquePartialR(1:4,4:6);
                        constraintTorquePJRI = constTorquePartialR(5:8,1:3);
                        constraintTorquePJRJ = constTorquePartialR(5:8,4:6);
                        
                        constraintTorquePIPI = constTorquePartialP(1:4,1:4);
                        constraintTorquePIPJ = constTorquePartialP(1:4,5:8);
                        constraintTorquePJPI = constTorquePartialP(5:8,1:4);
                        constraintTorquePJPJ = constTorquePartialP(5:8,5:8);
                        
                        % Place these partial derivatives in the correct
                        % location of the Jacobian
                        % Partial derivative of forces w.r.t. r
                        constForcesTotalPartialR((3*bodyInew - 2):(3*bodyInew),(3*bodyInew - 2):(3*bodyInew)) = ...
                            constForcesTotalPartialR((3*bodyInew - 2):(3*bodyInew),(3*bodyInew - 2):(3*bodyInew)) + constraintForceRIRI;
                        constForcesTotalPartialR((3*bodyInew - 2):(3*bodyInew),(3*bodyJnew - 2):(3*bodyJnew)) = ...
                            constForcesTotalPartialR((3*bodyInew - 2):(3*bodyInew),(3*bodyJnew - 2):(3*bodyJnew)) + constraintForceRIRJ;
                        constForcesTotalPartialR((3*bodyJnew - 2):(3*bodyJnew),(3*bodyInew - 2):(3*bodyInew)) = ...
                            constForcesTotalPartialR((3*bodyJnew - 2):(3*bodyJnew),(3*bodyInew - 2):(3*bodyInew)) + constraintForceRJRI;
                        constForcesTotalPartialR((3*bodyJnew - 2):(3*bodyJnew),(3*bodyJnew - 2):(3*bodyJnew)) = ...
                            constForcesTotalPartialR((3*bodyJnew - 2):(3*bodyJnew),(3*bodyJnew - 2):(3*bodyJnew)) + constraintForceRJRJ;
                        
                        % Partial derivative of forces w.r.t. p
                        constForcesTotalPartialP((3*bodyInew - 2):(3*bodyInew),(4*bodyInew - 3):(4*bodyInew)) =  ...
                            constForcesTotalPartialP((3*bodyInew - 2):(3*bodyInew),(4*bodyInew - 3):(4*bodyInew)) + constraintForceRIPI;
                        constForcesTotalPartialP((3*bodyInew - 2):(3*bodyInew),(4*bodyJnew - 3):(4*bodyJnew)) = ...
                            constForcesTotalPartialP((3*bodyInew - 2):(3*bodyInew),(4*bodyJnew - 3):(4*bodyJnew)) + constraintForceRIPJ;
                        constForcesTotalPartialP((3*bodyJnew - 2):(3*bodyJnew),(4*bodyInew - 3):(4*bodyInew)) = ...
                            constForcesTotalPartialP((3*bodyJnew - 2):(3*bodyJnew),(4*bodyInew - 3):(4*bodyInew)) + constraintForceRJPI;
                        constForcesTotalPartialP((3*bodyJnew - 2):(3*bodyJnew),(4*bodyJnew - 3):(4*bodyJnew)) = ...
                            constForcesTotalPartialP((3*bodyJnew - 2):(3*bodyJnew),(4*bodyJnew - 3):(4*bodyJnew)) + constraintForceRJPJ;
                        
                        % Partial derivative of torques w.r.t. r
                        constTorquesTotalPartialR((4*bodyInew - 3):(4*bodyInew),(3*bodyInew - 2):(3*bodyInew)) =  ...
                            constTorquesTotalPartialR((4*bodyInew - 3):(4*bodyInew),(3*bodyInew - 2):(3*bodyInew)) + constraintTorquePIRI;
                        constTorquesTotalPartialR((4*bodyInew - 3):(4*bodyInew),(3*bodyJnew - 2):(3*bodyJnew)) = ...
                            constTorquesTotalPartialR((4*bodyInew - 3):(4*bodyInew),(3*bodyJnew - 2):(3*bodyJnew)) + constraintTorquePIRJ;
                        constTorquesTotalPartialR((4*bodyJnew - 3):(4*bodyJnew),(3*bodyInew - 2):(3*bodyInew)) = ...
                            constTorquesTotalPartialR((4*bodyJnew - 3):(4*bodyJnew),(3*bodyInew - 2):(3*bodyInew)) + constraintTorquePJRI;
                        constTorquesTotalPartialR((4*bodyJnew - 3):(4*bodyJnew),(3*bodyJnew - 2):(3*bodyJnew)) = ...
                            constTorquesTotalPartialR((4*bodyJnew - 3):(4*bodyJnew),(3*bodyJnew - 2):(3*bodyJnew)) + constraintTorquePJRJ;
                        
                        % Partial derivative of torques w.r.t. p
                        constTorquesTotalPartialP((4*bodyInew - 3):(4*bodyInew),(4*bodyInew - 3):(4*bodyInew)) =  ...
                            constTorquesTotalPartialP((4*bodyInew - 3):(4*bodyInew),(4*bodyInew - 3):(4*bodyInew)) + constraintTorquePIPI;
                        constTorquesTotalPartialP((4*bodyInew - 3):(4*bodyInew),(4*bodyJnew - 3):(4*bodyJnew)) =  ...
                            constTorquesTotalPartialP((4*bodyInew - 3):(4*bodyInew),(4*bodyJnew - 3):(4*bodyJnew)) + constraintTorquePIPJ;
                        constTorquesTotalPartialP((4*bodyJnew - 3):(4*bodyJnew),(4*bodyInew - 3):(4*bodyInew)) = ...
                            constTorquesTotalPartialP((4*bodyJnew - 3):(4*bodyJnew),(4*bodyInew - 3):(4*bodyInew)) + constraintTorquePJPI;
                        constTorquesTotalPartialP((4*bodyJnew - 3):(4*bodyJnew),(4*bodyJnew - 3):(4*bodyJnew)) = ...
                            constTorquesTotalPartialP((4*bodyJnew - 3):(4*bodyJnew),(4*bodyJnew - 3):(4*bodyJnew)) + constraintTorquePJPJ;
                        
                    end
                else
                    % Include both bodyI and bodyJ in this case
                    constraintForceRIRI = constForcePartialR(1:3,1:3);
                    constraintForceRIRJ = constForcePartialR(1:3,4:6);
                    constraintForceRJRI = constForcePartialR(4:6,1:3);
                    constraintForceRJRJ = constForcePartialR(4:6,4:6);
                    
                    constraintForceRIPI = constForcePartialP(1:3,1:4);
                    constraintForceRIPJ = constForcePartialP(1:3,5:8);
                    constraintForceRJPI = constForcePartialP(4:6,1:4);
                    constraintForceRJPJ = constForcePartialP(4:6,5:8);
                    
                    constraintTorquePIRI = constTorquePartialR(1:4,1:3);
                    constraintTorquePIRJ = constTorquePartialR(1:4,4:6);
                    constraintTorquePJRI = constTorquePartialR(5:8,1:3);
                    constraintTorquePJRJ = constTorquePartialR(5:8,4:6);
                    
                    constraintTorquePIPI = constTorquePartialP(1:4,1:4);
                    constraintTorquePIPJ = constTorquePartialP(1:4,5:8);
                    constraintTorquePJPI = constTorquePartialP(5:8,1:4);
                    constraintTorquePJPJ = constTorquePartialP(5:8,5:8);
                    
                    % Place these partial derivatives in the correct
                    % location of the Jacobian
                    % Partial derivative of forces w.r.t. r
                    constForcesTotalPartialR((3*bodyI - 2):(3*bodyI),(3*bodyI - 2):(3*bodyI)) = ...
                        constForcesTotalPartialR((3*bodyI - 2):(3*bodyI),(3*bodyI - 2):(3*bodyI)) + constraintForceRIRI;
                    constForcesTotalPartialR((3*bodyI - 2):(3*bodyI),(3*bodyJ - 2):(3*bodyJ)) = ...
                        constForcesTotalPartialR((3*bodyI - 2):(3*bodyI),(3*bodyJ - 2):(3*bodyJ)) + constraintForceRIRJ;
                    constForcesTotalPartialR((3*bodyJ - 2):(3*bodyJ),(3*bodyI - 2):(3*bodyI)) = ...
                        constForcesTotalPartialR((3*bodyJ - 2):(3*bodyJ),(3*bodyI - 2):(3*bodyI)) + constraintForceRJRI;
                    constForcesTotalPartialR((3*bodyJ - 2):(3*bodyJ),(3*bodyJ - 2):(3*bodyJ)) = ...
                        constForcesTotalPartialR((3*bodyJ - 2):(3*bodyJ),(3*bodyJ - 2):(3*bodyJ)) + constraintForceRJRJ;
                    
                    % Partial derivative of forces w.r.t. p
                    constForcesTotalPartialP((3*bodyI - 2):(3*bodyI),(4*bodyI - 3):(4*bodyI)) =  ...
                        constForcesTotalPartialP((3*bodyI - 2):(3*bodyI),(4*bodyI - 3):(4*bodyI)) + constraintForceRIPI;
                    constForcesTotalPartialP((3*bodyI - 2):(3*bodyI),(4*bodyJ - 3):(4*bodyJ)) = ...
                        constForcesTotalPartialP((3*bodyI - 2):(3*bodyI),(4*bodyJ - 3):(4*bodyJ)) + constraintForceRIPJ;
                    constForcesTotalPartialP((3*bodyJ - 2):(3*bodyJ),(4*bodyI - 3):(4*bodyI)) = ...
                        constForcesTotalPartialP((3*bodyJ - 2):(3*bodyJ),(4*bodyI - 3):(4*bodyI)) + constraintForceRJPI;
                    constForcesTotalPartialP((3*bodyJ - 2):(3*bodyJ),(4*bodyJ - 3):(4*bodyJ)) = ...
                        constForcesTotalPartialP((3*bodyJ - 2):(3*bodyJ),(4*bodyJ - 3):(4*bodyJ)) + constraintForceRJPJ;
                    
                    % Partial derivative of torques w.r.t. r
                    constTorquesTotalPartialR((4*bodyI - 3):(4*bodyI),(3*bodyI - 2):(3*bodyI)) =  ...
                        constTorquesTotalPartialR((4*bodyI - 3):(4*bodyI),(3*bodyI - 2):(3*bodyI)) + constraintTorquePIRI;
                    constTorquesTotalPartialR((4*bodyI - 3):(4*bodyI),(3*bodyJ - 2):(3*bodyJ)) = ...
                        constTorquesTotalPartialR((4*bodyI - 3):(4*bodyI),(3*bodyJ - 2):(3*bodyJ)) + constraintTorquePIRJ;
                    constTorquesTotalPartialR((4*bodyJ - 3):(4*bodyJ),(3*bodyI - 2):(3*bodyI)) = ...
                        constTorquesTotalPartialR((4*bodyJ - 3):(4*bodyJ),(3*bodyI - 2):(3*bodyI)) + constraintTorquePJRI;
                    constTorquesTotalPartialR((4*bodyJ - 3):(4*bodyJ),(3*bodyJ - 2):(3*bodyJ)) = ...
                        constTorquesTotalPartialR((4*bodyJ - 3):(4*bodyJ),(3*bodyJ - 2):(3*bodyJ)) + constraintTorquePJRJ;
                    
                    % Partial derivative of torques w.r.t. p
                    constTorquesTotalPartialP((4*bodyI - 3):(4*bodyI),(4*bodyI - 3):(4*bodyI)) =  ...
                        constTorquesTotalPartialP((4*bodyI - 3):(4*bodyI),(4*bodyI - 3):(4*bodyI)) + constraintTorquePIPI;
                    constTorquesTotalPartialP((4*bodyI - 3):(4*bodyI),(4*bodyJ - 3):(4*bodyJ)) =  ...
                        constTorquesTotalPartialP((4*bodyI - 3):(4*bodyI),(4*bodyJ - 3):(4*bodyJ)) + constraintTorquePIPJ;
                    constTorquesTotalPartialP((4*bodyJ - 3):(4*bodyJ),(4*bodyI - 3):(4*bodyI)) = ...
                        constTorquesTotalPartialP((4*bodyJ - 3):(4*bodyJ),(4*bodyI - 3):(4*bodyI)) + constraintTorquePJPI;
                    constTorquesTotalPartialP((4*bodyJ - 3):(4*bodyJ),(4*bodyJ - 3):(4*bodyJ)) = ...
                        constTorquesTotalPartialP((4*bodyJ - 3):(4*bodyJ),(4*bodyJ - 3):(4*bodyJ)) + constraintTorquePJPJ;
                    
                end
            end
            obj.myConstraintForcePartialR = constForcePartialR;
            obj.myConstraintForcePartialP = constForcePartialP;
            obj.myConstraintTorquePartialR = constTorquePartialR;
            obj.myConstraintTorquePartialP = constTorquePartialP;
        end
        function obj = computeQuasiNewtonPsi(obj)
            % Compute the iteration matrix (psi) for the quasi-newton
            % method
            
            % Extract terms for iteration matrix. All of the terms were
            % previously compute when we computed the Gmatrix, so they do
            % not need to be recomputed.
            massMatrix = obj.myMassMatrix;
            JpMatrix = obj.myJpMatrixTotal;
            Pmatrix = obj.myPMatrix;
            phiPartialR = obj.myPhiPartialR;
            phiPartialP = obj.myPhiPartialP;
            nB = obj.myNumBodiesMinusGround;
            nC = obj.myNumConstraints;
            
            % Populate components of psi matrix
            fatM = [massMatrix, zeros(3*nB,4*nB), zeros(3*nB,nB);
                zeros(4*nB,3*nB), JpMatrix, Pmatrix';
                zeros(nB,3*nB), Pmatrix, zeros(nB,nB)];
            fatC = [phiPartialR, phiPartialP, zeros(nC,nB)];
            
            % Matrix Assemble!!!
            psi = [fatM fatC';
                fatC zeros(nC,nC)];
            obj.myPsi = psi;
            
        end
        function obj = computeGMatrix(obj, stepsize, order)
            % Compute the g-matrix for the BDF method iteration
            %
            % Function inputs:
            % stepsize : double
            %   Stepsize used for numerical integration.
            %
            % order : int
            %   Order that was used for the current step of the BDF method
            
            % Set beta0 depending on the order of the BDF method being used.
            if (order == 1)
                beta0 = 1;
            elseif (order == 2)
                beta0 = 2/3;
            end
            
            % Compute and extract the properties needed for the g matrix
            massMatrix = obj.myMassMatrix;
            
            constLambdaGuess = obj.myConstraintLagrangeMultipliers;
            pLambdaGuess = obj.myEulerParamLagrangeMultipliers;
            
            if (obj.myBodyIsGround == 1)
                rDDotGuess = obj.myRDDot(4:end,1);
                pDDotGuess = obj.myPDDot(5:end,1);
            else
                rDDotGuess = obj.myRDDot;
                pDDotGuess = obj.myPDDot;
            end
            
            obj.computePhiPartialR();
            phiPartialR = obj.myPhiPartialR;
            
            obj.computePhiPartialP();
            phiPartialP = obj.myPhiPartialP;
            
            obj.computeForceVector();
            forceVector = obj.myForceVector;
            
            obj.computeJpMatrixTotal();
            JpMatrix = obj.myJpMatrixTotal;
            
            obj.computePMatrix();
            Pmatrix = obj.myPMatrix;
            
            obj.computeTorqueHatVector();
            tauHat = obj.myTorqueHatVector;
            
            obj.computePhiP();
            phiP = obj.myPhiP;
            
            obj.computePhiK();
            phiK = obj.myPhiK;
            
            obj.computePhiD();
            phiD = obj.myPhiD;
            
            phi = [phiK;
                phiD];
            
            % Compute all terms of Gmatrix
            forceTerm = massMatrix*rDDotGuess + phiPartialR'*constLambdaGuess - forceVector;
            torqueTerm = JpMatrix*pDDotGuess + phiPartialP'*constLambdaGuess + Pmatrix'*pLambdaGuess - tauHat;
            eulerParamTerm = 1/(beta0^2*stepsize^2)*phiP;
            constraintTerm = 1/(beta0^2*stepsize^2)*phi;
            %             eulerParamTerm = phiP;
            %             constraintTerm = phi;
            
            
            % Matrix Assemble!!!
            gMatrix = [forceTerm;
                torqueTerm;
                eulerParamTerm;
                constraintTerm];
            
            obj.myGMatrix = gMatrix;
        end
        function obj = BDFmethodIteration(obj, rDDotGuess, pDDotGuess, order, stepsize, time, stepNumber)
            % Perform one iteration of the BDF method
            %
            % Function inputs:
            % rDDotGuess : 3N x 1 vec, where N is the number of bodies not including the ground.
            %   Vector containing the guess for the position accelerations of each
            %   body in the system, not including the ground.
            %
            % pDDotGuess : 4N x 1 vec, where N is the number of bodies not including the ground.
            %   Vector containing the guess for the second time derivative
            %   of the Euler parameters of each body in the system, not
            %   including the ground.
            %
            % order : int
            %   Order of the BDF method. Currently only support 1 and 2
            %
            % stepsize : double
            %   Stepsize of the numerical integration
            %
            % time : double
            %   Current time of the system
            %
            % stepNumber : int
            %   Current step number of the dynamics analysis.
            
            if (order == 1)
                % Extract the solution from the previous time step. If
                % there is a ground, remove it from the vector.
                if (obj.myBodyIsGround == 1)
                    rNminus1 = obj.myRTotal(4:end,(stepNumber - 1));
                    rDotNminus1 = obj.myRDotTotal(4:end,(stepNumber - 1));
                    
                    pNminus1 = obj.myPTotal(5:end,(stepNumber - 1));
                    pDotNminus1 = obj.myPDotTotal(5:end,(stepNumber - 1));
                else
                    rNminus1 = obj.myRTotal(:,(stepNumber - 1));
                    rDotNminus1 = obj.myRDotTotal(:,(stepNumber - 1));
                    
                    pNminus1 = obj.myPTotal(:,(stepNumber - 1));
                    pDotNminus1 = obj.myPDotTotal(:,(stepNumber - 1));
                end
                
                % Compute next guess for rDot and pDot.
                rDotGuess = rDotNminus1 + stepsize*rDDotGuess;
                pDotGuess = pDotNminus1 + stepsize*pDDotGuess;
                
                % Compute next guess for r and p.
                rGuess = rNminus1 + stepsize*rDotGuess;
                pGuess = pNminus1 + stepsize*pDotGuess;
                
            elseif (order == 2)
                % Extract the solution from the previous two time steps.
                % Remove ground if it is present in the system.
                if (obj.myBodyIsGround == 1)
                    rNminus1 = obj.myRTotal(4:end,(stepNumber - 1));
                    rNminus2 = obj.myRTotal(4:end,(stepNumber - 2));
                    rDotNminus1 = obj.myRDotTotal(4:end,(stepNumber - 1));
                    rDotNminus2 = obj.myRDotTotal(4:end,(stepNumber - 2));
                    
                    pNminus1 = obj.myPTotal(5:end,(stepNumber - 1));
                    pNminus2 = obj.myPTotal(5:end,(stepNumber - 2));
                    pDotNminus1 = obj.myPDotTotal(5:end,(stepNumber - 1));
                    pDotNminus2 = obj.myPDotTotal(5:end,(stepNumber - 2));
                else
                    rNminus1 = obj.myRTotal(:,(stepNumber - 1));
                    rNminus2 = obj.myRTotal(:,(stepNumber - 2));
                    rDotNminus1 = obj.myRDotTotal(:,(stepNumber - 1));
                    rDotNminus2 = obj.myRDotTotal(:,(stepNumber - 2));
                    
                    pNminus1 = obj.myPTotal(:,(stepNumber - 1));
                    pNminus2 = obj.myPTotal(:,(stepNumber - 2));
                    pDotNminus1 = obj.myPDotTotal(:,(stepNumber - 1));
                    pDotNminus2 = obj.myPDotTotal(:,(stepNumber - 2));
                end
                % Compute next guess for rDot and pDot.
                rDotGuess = (4/3)*rDotNminus1 - (1/3)*rDotNminus2 + (2/3)*stepsize*rDDotGuess;
                pDotGuess = (4/3)*pDotNminus1 - (1/3)*pDotNminus2 + (2/3)*stepsize*pDDotGuess;
                
                % Compute next guess for r and p.
                rGuess = (4/3)*rNminus1 - (1/3)*rNminus2 + (2/3)*stepsize*rDotGuess;
                pGuess = (4/3)*pNminus1 - (1/3)*pNminus2 + (2/3)*stepsize*pDotGuess;
                
            else
                error('Only BDF method of order 1 or 2 supported')
            end
            
            % If there is a ground in the system, add it back in.
            if (obj.myBodyIsGround == 1)
                rGuess = [zeros(3,1); rGuess];
                rDotGuess = [zeros(3,1); rDotGuess];
                rDDotGuess = [zeros(3,1); rDDotGuess];
                
                pGuess = [zeros(4,1); pGuess];
                pDotGuess = [zeros(4,1); pDotGuess];
                pDDotGuess = [zeros(4,1); pDDotGuess];
            end
            
            % Store the guesses for the positions, velocities, and
            % accelerations as the current state of the system.
            obj.updateSystemState(rGuess, rDotGuess, rDDotGuess, pGuess, pDotGuess, pDDotGuess, time);
        end
        function obj = solveFirstStepDynamicsAnalysis(obj, time)
            % Compute the second time derivative of the generalized
            % coordinates and the Lagrange multipliers for the constraint
            % forces for the first time step in a dynamics analysis.
            
            % Compute left hand side of equations of motion
            obj.computeLHSforEOM();
            LHSforEOM = obj.myLHSforEOM;
            
            % Compute right hand side of equations of motion
            obj.computeRHSforEOM();
            RHSforEOM = obj.myRHSforEOM;
            
            % Solve for accelerations
            solution = LHSforEOM\RHSforEOM;
            
            % Divide the solution matrix into qDDot and the Lagrange
            % multipliers
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            nConst = obj.myNumConstraints;
            
            qDDot = solution(1:7*nBodies);
            lagrangeMultsP = solution((7*nBodies + 1):(8*nBodies));
            constraintLagrangeMults = solution((8*nBodies + 1):(8*nBodies + nConst));
            
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
            
            % Update the accelerations
            obj.updateSystemState([], [], rDDot, [], [], pDDot, time);
            
            % Store the lagrange multipliers
            obj.myConstraintLagrangeMultipliers = constraintLagrangeMults;
            obj.myEulerParamLagrangeMultipliers = lagrangeMultsP;
        end
        function obj = computeLHSforEOM(obj)
            % Compute the left hand side for the matrix representation of the
            % Newton-Euler form of the equations of motion.
            
            % Extract mass matrix
            massMatrix = obj.myMassMatrix;
            
            % Compute Jp matrix
            obj.computeJpMatrixTotal();
            JpMatrix = obj.myJpMatrixTotal;
            
            % Extract Pmatrix, phiPartialR, and phiPartialP. These were
            % computed when performing the initial condition check in the
            % previous step so save some time by not recomputing them.
            phiPartialR = obj.myPhiPartialR;
            phiPartialP = obj.myPhiPartialP;
            Pmatrix = obj.myPMatrix;
            
            % Determine number of bodies and number of constraints
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            nConst = obj.myNumConstraints;
            
            % Matrix Assemble!!!
            LHSmatrix = [massMatrix, zeros(3*nBodies,4*nBodies), zeros(3*nBodies,nBodies), phiPartialR';
                zeros(4*nBodies,3*nBodies), JpMatrix, Pmatrix', phiPartialP';
                zeros(nBodies,3*nBodies), Pmatrix, zeros(nBodies,nBodies), zeros(nBodies,nConst);
                phiPartialR, phiPartialP, zeros(nConst,nBodies), zeros(nConst,nConst)];
            obj.myLHSforEOM = LHSmatrix;
        end
        function obj = computeRHSforEOM(obj)
            % Compute the right hand side for the matrix-form of the
            % Newton-Euler equations of motion
            
            % Compute the total external force acting on the system
            obj.computeForceVector();
            forceVector = obj.myForceVector;
            
            % Compute total external torque acting on the system
            obj.computeTorqueHatVector();
            tauHat = obj.myTorqueHatVector;
            
            % Compute gamma total, which contains gamma for the Euler
            % parameter normalization constraints and gamma for the
            % kinematic and driving constraints.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            nConst = obj.myNumConstraints;
            obj.computeGamma();
            gamma = obj.myGamma;
            gammaHat = gamma(1:nConst);
            gammaP = gamma((nConst+1):(nConst+nBodies));
            
            % Assemble RHS vector
            RHSvector = [forceVector;
                tauHat;
                gammaP;
                gammaHat];
            obj.myRHSforEOM = RHSvector;
        end
        function obj = checkInitialConditions(obj)
            % Check initial conditions of a dynamics analysis to make sure
            % they satisfy the constraint equations
            
            % Compute phi full and check to make sure it equals (approximately) zero.
            obj.computePhiFull();
            phiK = obj.myPhiK;
            phiD = obj.myPhiD;
            phiKD = [phiK; phiD];
            if (abs(norm(phiKD)) > 10^-3)
                error('Constraint matrix not equal to zero in initial pose.');
                phiFullFlag = 0;
            else
                phiFullFlag = 1;
            end
            
            % Check Euler parameter normalization constraint to make sure
            % it is zero. PhiP was already computed as part of phiFull so
            % you can just pull it.
            if (abs(norm(obj.myPhiP)) > 10^-2)
                error('Euler parameter normalization constraint not satisified in initial pose.')
                eulerParamConstFlag = 0;
            else
                eulerParamConstFlag = 1;
            end
            
            % Check the velocity constraint.
            % Compute necessary parameters.
            obj.computePhiPartialR();
            obj.computePhiPartialP();
            obj.computeNu();
            
            % Extract necessary parameters.
            phiPartialR = obj.myPhiPartialR;
            phiPartialP = obj.myPhiPartialP;
            nu = obj.myNu;
            rDot = obj.myRDot;
            pDot = obj.myPDot;
            
            % Remove ground from system if it exists
            if (obj.myBodyIsGround == 1)
                % rDot and pDot include the ground so remove them.
                rDotVec = rDot(4:end,1);
                pDotVec = pDot(5:end,1);
            else
                rDotVec = rDot;
                pDotVec = pDot;
            end
            
            % Extract the portion of nu that relates to the constraints
            nuConst = nu(1:obj.myNumConstraints);
            
            % Check to ensure this constraint is satisfied
            check = phiPartialR*rDotVec + phiPartialP*pDotVec - nuConst;
            if (abs(norm(check)) > 10^-6)
                error('Velocity constraint not satisified in initial pose.')
                velocityConstFlag = 0;
            else
                velocityConstFlag = 1;
            end
            
            % Check that the Euler parameter velocity normalization
            % constraint is satisfied.
            obj.computePMatrix();
            Pmatrix = obj.myPMatrix;
            
            check2 = Pmatrix*pDotVec;
            if (abs(norm(check2)) > 10^-7)
                error('Euler parameter velocity normalization constraint is not satisified in initial pose.')
                eulerParamVelConstFlag = 0;
            else
                eulerParamVelConstFlag = 1;
            end
            
            % Check to see if any of the flags are zero. If they are, one
            % of the constraints was not satisfied
            if (phiFullFlag == 0) || (eulerParamConstFlag == 0) || (velocityConstFlag == 0) || (eulerParamVelConstFlag == 0)
                obj.myInitCondFlag = 0;
            else
                obj.myInitCondFlag = 1;
            end
        end
        function obj = assemblyAnalysis(obj)
            % Perform an assembly analysis before performing any other
            % sort of analysis (e.g. kinematics, inverse dynamics, or
            % dynamics)
            
            % Set parameters
            rVal = 1;
            maxIter = 50;
            opt = optimoptions('fminunc');
            opt.GradObj = 'on';
            tol = 10^-8;
            
            % Extract initial values
            if (obj.myBodyIsGround == 1)
                r0 = obj.myR(4:end);
                p0 = obj.myP(5:end);
            else
                r0 = obj.myR;
                p0 = obj.myP;
            end
            q0 = [r0; p0];
            
            
            iter = 1;
            while iter <= maxIter
                % Use fminunc to minimize psi
                [qMin, fMin, exitflag, output] = fminunc(@(q)obj.computeAssemblyAnalysisPsi(q,q0,rVal), q0, opt);
                
                % Check to see if this new value of q has converged
                deltaQ = qMin - q0;
                if (norm(deltaQ) < tol)
                    break;
                end
                
                % If q has not converged, increment r and continue the
                % process.
                q0 = qMin;
                rVal = rVal*2;

                iter = iter + 1;
            end
            
            % Update the system state with these final values of r and p.
            nBodies = obj.myNumBodiesMinusGround;
            if (obj.myBodyIsGround == 1)
                rTemp = qMin(1:3*nBodies);
                pTemp = qMin((3*nBodies+1):(7*nBodies));
                
                rFinal = [0; 0; 0; rTemp];
                pFinal = [1; 0; 0; 0; pTemp];
            else
                rFinal = qMin(1:3*nBodies);
                pFinal = qMin((3*nBodies+1):(7*nBodies));
            end
            
            t = 0; % Because assembly analysis is only performed at t = 0.
            obj.updateSystemState(rFinal, [], [], pFinal, [], [], t);            
        end
        function [psi, psiPartialQ] = computeAssemblyAnalysisPsi(obj, q, q0, r)
            % Compute the value of psi and the partial derivative of psi
            % w.r.t. q for the assembly analysis.
            %
            % Function inputs:
            % q : 7*nBodies x 1 vector
            %   Current guess for the value of q that minimizes psi.
            %
            % q0 : 7*nBodies x 1 vector
            %   Initial guess for the value of q.
            %
            % r : double
            %   Current value of r used in computing psi and psiPartialQ.
            
            % Set the system values to the current guess of q. 
            nBodies = obj.myNumBodiesMinusGround;
            if (obj.myBodyIsGround == 1)
                rTemp = q(1:3*nBodies);
                pTemp = q((3*nBodies+1):(7*nBodies));
                
                rGuess = [0; 0; 0; rTemp];
                pGuess = [1; 0; 0; 0; pTemp];
            else
                rGuess = q(1:3*nBodies);
                pGuess = q((3*nBodies+1):(7*nBodies));
            end
            
            t = 0; % Assembly analysis is only run at t = 0;
            obj.updateSystemState(rGuess, [], [], pGuess, [], [], t);
            
            % Compute the components of psi and psiPartialQ
            obj.computePhiFull();
            phiFull = obj.myPhiFull();
            
            obj.computePhiFullJacobian();
            phiFullJacobian = obj.myPhiFullJacobian();
            
            % Compute psi and psiPartialQ.
            psi = (q - q0)'*(q - q0) + r*(phiFull'*phiFull);
            psiPartialQ = 2*(q - q0)' + 2*r*phiFull'*phiFullJacobian;           
        end
        function obj = inverseDynamicsAnalysis(obj, startTime, endTime, timestep, displayFlag)
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
            %
            % displayFlag : int
            %   Flag for user to indicate if they want to display when each
            %   time step of analysis has been completed.
            
            % Check to make sure this is a fully constrained system
            nConstTotal = obj.myNumConstraints + obj.myNumBodiesMinusGround;
            if (nConstTotal < 7*obj.myNumBodiesMinusGround)
                error('System not fully constrained. Cannot perform inverse dynamics analysis')
            elseif (nConstTotal > 7*obj.myNumBodiesMinusGround)
                error('System over-constrained. Cannot perform inverse dynamics analysis')
            end
            
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
                obj.storeSystemState(iT);
                
                % Compute the Lagrange multipliers using the equations of
                % motion
                obj.computeLagrangeMultipliers();
                
                % Compute forces and torques associated with each Lagrange
                % multiplier.
                obj.computeConstraintForces();
                obj.computeConstraintTorques();
                
                % Store the current forces and torques for all the
                % constraints
                obj.storeConstraintForcesAndTorques(t,iT);
                
                if (displayFlag == 1)
                    disp(['Inverse dynamics analysis completed for ' num2str(t) ' sec']);
                end
            end
            
            
        end
        function obj = storeConstraintForcesAndTorques(obj,t,iT)
            % Store forces and torques that were computed for each constraint
            % at this specific time step. The format of both of this
            % matrices is 3nc x nTimesteps, where every three rows
            % represents the forces or torques for a different constraint
            % and each column is a different timestep.
            nTimeSteps = obj.myNumTimeSteps;
            nConst = obj.myNumConstraints;
            nBodies = obj.myNumBodies;
            
            % Store time
            obj.myTimeTotal(1,iT) = t;
            
            % For each body, store the force and torque for this specific
            % time step.
            for iB = 1:nBodies
                force = obj.myBodies{iB}.myConstraintForces;
                torque = obj.myBodies{iB}.myConstraintTorques;
                torqueOmega = obj.myBodies{iB}.myConstraintTorquesOmega;
                forceVec = reshape(force,[3*nConst,1]);
                torqueVec = reshape(torque,[4*nConst,1]);
                torqueOmegaVec = reshape(torqueOmega,[3*nConst,1]);
                
                % Store forces and torques
                obj.myBodies{iB}.myConstraintForcesTotal(:,iT) = forceVec;
                obj.myBodies{iB}.myConstraintTorquesTotal(:,iT) = torqueVec;
                obj.myBodies{iB}.myConstraintTorquesOmegaTotal(:,iT) = torqueOmegaVec;
            end
        end
        function obj = computeConstraintForces(obj)
            % Compute the forces induced by each constraint.
            nConst = obj.myNumConstraints;
            nBodies = obj.myNumBodies;
            lagrangeMults = obj.myConstraintLagrangeMultipliers;
            
            % Need to compute phiPartialR for the final values of R
            obj.computePhiPartialR();
            phiPartialR = obj.myPhiPartialR;
            
            % Seed the constraint forces for all bodies
            for iB = 1:nBodies
                obj.myBodies{iB}.myConstraintForces = zeros(3,nConst);
            end
            
            % For all bodies, extract the correct portion of
            % phiPartialR, compute the force due to constraint iC on
            % this body, and store this force at the body level. Ignore
            % the body if it is the ground. The constraint torque will
            % just remain all zeros for the ground.
            for iC = 1:nConst
                count = 1;
                for iB = 1:nBodies
                    if (obj.myBodies{iB}.myIsGround == 0)
                        % Compute the constraint forces
                        constraintForce = -phiPartialR(iC,(3*count-2):3*count)'*lagrangeMults(iC);
                        obj.myBodies{iB}.myConstraintForces(:,iC) = constraintForce;
                        count = count + 1;
                    end
                end
            end
        end
        function obj = computeConstraintTorques(obj)
            % Compute the torques induced by each constraint
            nConst = obj.myNumConstraints;
            nBodies = obj.myNumBodies;
            lagrangeMults = obj.myConstraintLagrangeMultipliers;
            
            % Compute phiPartial P for the final values of P.
            obj.computePhiPartialP();
            phiPartialP = obj.myPhiPartialP;
            
            % Seed the constraint forces for all bodies
            for iB = 1:nBodies
                obj.myBodies{iB}.myConstraintTorques = zeros(4,nConst);
            end
            
            % For all bodies, extract the correct portion of
            % phiPartialP, compute the torque due to constraint iC on
            % this body, and store this torque at the body level. Ignore
            % the body if it is the ground. The constraint torque will
            % just remain all zeros for the ground.
            for iC = 1:nConst
                count = 1;
                for iB = 1:nBodies
                    if (obj.myBodies{iB}.myIsGround == 0)
                        % Compute the constraint forces
                        constraintTorque = -phiPartialP(iC,(4*count-3):4*count)'*lagrangeMults(iC);
                        obj.myBodies{iB}.myConstraintTorques(:,iC) = constraintTorque;
                        count = count + 1;
                    end
                end
            end
            
            % Convert constraint torques from r-p formulation to r-omega
            % formulation for each body
            for iB = 1:nBodies
                obj.myBodies{iB}.convertConstraintTorques();
            end
        end
        function obj = computeLagrangeMultipliers(obj)
            % Compute the Lagrange multipliers that will be used for
            % computing the reaction force and torques
            
            % Compute the RHS for linear system of equations for inverse dynamics analysis
            obj.computeInvDynRHS();
            invDynRHS = obj.myInvDynRHS;
            
            % Compute the constraint partial derivative matrix. The matrix
            % is the A matrix in solving Ax = b.
            obj.computeInvDynMatrix();
            invDynMatrix = obj.myInvDynMatrix();
            
            
            % Compute Lagrange multipliers
            lagrangeMultipliers = invDynMatrix\invDynRHS;
            
            % Extract lagrange multipliers for constraints
            nConst = obj.myNumConstraints;
            lagrangeMultConst = lagrangeMultipliers(1:nConst);
            lagrangeMultEulerParams = lagrangeMultipliers((nConst+1):end);
            obj.myConstraintLagrangeMultipliers = lagrangeMultConst;
            obj.myEulerParamLagrangeMultipliers = lagrangeMultEulerParams;
            
        end
        function obj = computeInvDynMatrix(obj)
            % Compute the left hand side matrix for the inverse dynamics
            % analysis. This matrix is of the form:
            % [phiPartialR' zeros(3*nb,nb);
            % phiPartialP' P']
            
            % Compute partial derivative of the constraint matrix w.r.t.
            % displacements.
            obj.computePhiPartialR();
            phiPartialR = obj.myPhiPartialR;
            
            % Compute partial derivative of the constraint matrix w.r.t.
            % orientation.
            obj.computePhiPartialP();
            phiPartialP = obj.myPhiPartialP;
            
            % Compute P matrix. This matrix contains the Euler
            % normalization constraints.
            obj.computePMatrix();
            P = obj.myPMatrix;
            
            % Build the A matrix
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            zeroMatrix = zeros(3*nBodies,nBodies);
            invDynMatrix = [phiPartialR', zeroMatrix;
                phiPartialP', P'];
            obj.myInvDynMatrix = invDynMatrix;
        end
        function obj = computePMatrix(obj)
            % Compute Euler param matrix
            
            % Determine number of bodies. Also, extract the Euler params
            % for all the bodies. Remove the first column of the parameters
            % if one of the bodies is that ground because we do not want
            % those values.
            p = obj.myP;
            pMatrixForm = reshape(p, [4,obj.myNumBodies]);
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
                pMatrixForm(:,1) = [];
            else
                nBodies = obj.myNumBodies;
            end
            
            % Seed the P matrix
            Pmatrix = zeros(nBodies,4*nBodies);
            
            
            % Loop through each body and add it to the matrix
            for iB = 1:nBodies
                Pmatrix(iB, (4*iB - 3):4*iB) = pMatrixForm(:,iB)';
            end
            obj.myPMatrix = Pmatrix;
        end
        function obj = computePhiPartialR(obj)
            % Compute partial derivative of the constraint matrix w.r.t
            % translation (R)
            
            % Seed the Jacobian. There will be 3 less columns in the
            % Jacobian if one of the bodies is the ground.
            nBodies = obj.myNumBodies;
            nKDconst = obj.myNumConstraints;
            
            if (obj.myBodyIsGround == 1)
                phiPartialR = zeros(nKDconst, 3*(nBodies - 1));
            else
                phiPartialR = zeros(nKDconst, 3*nBodies);
            end
            
            % Loop through each kinematic and driving constraint. Compute
            % phiPartialR for each. Insert these the
            % correct location in the Jacobian. This location depends on
            % which bodies are involved in the constraint. If one of the
            % bodies is the ground, then it will not appear in the
            % Jacobian.
            for iC = 1:nKDconst
                % Extract the bodies for this constraint
                bodyI = obj.myConstraints{iC}.myBodyI;
                bodyJ = obj.myConstraints{iC}.myBodyJ;
                
                % Compute phiPartialR
                time = obj.myTime;
                phiPartialRFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, 0, phiPartialRFlag, 0, 0, 0, 0, 0);
                phiR = obj.myConstraints{iC}.myPhiPartialR;
                
                % If one of the bodies in the system is the ground adjust
                % the body numbers accordingly.
                if (obj.myBodyIsGround == 1)
                    if (obj.myBodies{bodyI}.myIsGround == 1)
                        % Only care about bodyJ in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. This only works because I
                        % am forcing the ground to be body 1.
                        phiPartialR(iC,(3*(bodyJ-1) - 2):3*(bodyJ-1)) = phiR;
                        
                    elseif (obj.myBodies{bodyJ}.myIsGround == 1)
                        % Only care about bodyI in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian.
                        phiPartialR(iC,(3*(bodyI-1) - 2):3*(bodyI-1)) = phiR;
                        
                    else
                        % Include both bodyI and bodyJ in this case because
                        % neither is the ground, but still decrease body
                        % number by 1.
                        % Populate the partial derivative w.r.t. position.
                        phiPartialR(iC,(3*(bodyI-1) - 2):3*(bodyI-1)) = phiR(1:3);
                        phiPartialR(iC,(3*(bodyJ-1) - 2):3*(bodyJ-1)) = phiR(4:6);
                    end
                else
                    % Populate the partial derivative w.r.t. position.
                    phiPartialR(iC,(3*bodyI - 2):3*bodyI) = phiR(1:3);
                    phiPartialR(iC,(3*bodyJ - 2):3*bodyJ) = phiR(4:6);
                end
            end
            obj.myPhiPartialR = phiPartialR;
            
        end
        function obj = computePhiPartialP(obj)
            % Compute partial derivative of the constraint matrix w.r.t
            % orientation (P)
            
            % Seed the Jacobian. There will be 4 less columns in the
            % Jacobian if one of the bodies is the ground.
            nBodies = obj.myNumBodies;
            nKDconst = obj.myNumConstraints;
            
            if (obj.myBodyIsGround == 1)
                phiPartialP = zeros(nKDconst, 4*(nBodies - 1));
            else
                phiPartialP = zeros(nKDconst, 4*nBodies);
            end
            
            % Loop through each kinematic and driving constraint. Compute
            % phiPartialR for each. Insert these the
            % correct location in the Jacobian. This location depends on
            % which bodies are involved in the constraint. If one of the
            % bodies is the ground, then it will not appear in the
            % Jacobian.
            for iC = 1:nKDconst
                % Extract the bodies for this constraint
                bodyI = obj.myConstraints{iC}.myBodyI;
                bodyJ = obj.myConstraints{iC}.myBodyJ;
                
                % Compute phiPartialR
                time = obj.myTime;
                phiPartialPFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, 0, 0, phiPartialPFlag, 0, 0, 0, 0);
                phiP = obj.myConstraints{iC}.myPhiPartialP;
                
                % If one of the bodies in the system is the ground adjust
                % the body numbers accordingly.
                if (obj.myBodyIsGround == 1)
                    if (obj.myBodies{bodyI}.myIsGround == 1)
                        % Only care about bodyJ in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. This only works because I
                        % am forcing the ground to be body 1.
                        phiPartialP(iC,(4*(bodyJ-1) - 3):4*(bodyJ-1)) = phiP;
                        
                    elseif (obj.myBodies{bodyJ}.myIsGround == 1)
                        % Only care about bodyI in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian.
                        phiPartialP(iC,(4*(bodyI-1) - 3):4*(bodyI-1)) = phiP;
                        
                    else
                        % Include both bodyI and bodyJ in this case because
                        % neither is the ground, but still decrease body
                        % number by 1.
                        % Populate the partial derivative w.r.t. position.
                        phiPartialP(iC,(4*(bodyI-1) - 3):4*(bodyI-1)) = phiP(1:4);
                        phiPartialP(iC,(4*(bodyJ-1) - 3):4*(bodyJ-1)) = phiP(5:8);
                    end
                else
                    % Populate the partial derivative w.r.t. position.
                    phiPartialP(iC,(4*bodyI - 3):4*bodyI) = phiP(1:4);
                    phiPartialP(iC,(4*bodyJ - 3):4*bodyJ) = phiP(5:8);
                end
            end
            obj.myPhiPartialP  = phiPartialP;
        end
        
        function obj = computeInvDynRHS(obj)
            % Compute the RHS of the linear system of equations for the
            % inverse dynamics analysis
            
            % Determine number of bodies, not including ground.
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            
            % Obtain force vector
            obj.computeForceVector();
            forceVec = obj.myForceVector;
            
            % Obtain mass matrix
            massMatrix = obj.myMassMatrix;
            
            % Obtain vector of acceleration and second time derivative of
            % Euler params for all bodies
            rDDot = obj.myRDDot;
            pDDot = obj.myPDDot;
            if (obj.myBodyIsGround == 1)
                % Remove ground from rDDot and pDDot matrices, if there is a ground
                % defined in the system. The ground is always body 1.
                rDDotVec = rDDot(4:end,1);
                pDDotVec = pDDot(5:end,1);
            end
            
            % Obtain torqueHat vector
            obj.computeTorqueHatVector();
            torqueHatVec = obj.myTorqueHatVector;
            
            % Obtain JpMatrix
            obj.computeJpMatrixTotal();
            JpMatrix = obj.myJpMatrixTotal;
            
            % Compute force related terms of RHS
            forceTerms = massMatrix*rDDotVec - forceVec;
            
            % Compute torque related terms of RHS
            torqueTerms = JpMatrix*pDDotVec - torqueHatVec;
            
            % Assemble RHS vector.
            invDynRHS = -1*[forceTerms;
                torqueTerms];
            
            obj.myInvDynRHS = invDynRHS;
        end
        function obj = computeTorqueHatVector(obj)
            % Compute torqueHat for all bodies in the system and assemble
            % into a vector.
            
            % Seed torqueHat vector
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            torqueHatVector = zeros(4*nBodies,1);
            
            % Loop through all bodies. Compute torqueHat. Add this
            % torqueHat to the overall vector
            iT = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeTorqueHat();
                    torqueHat = obj.myBodies{iB}.myTorqueHat;
                    torqueHatVector((4*iT-3):(4*iT),1) = torqueHat;
                    iT = iT + 1;
                end
            end
            obj.myTorqueHatVector = torqueHatVector;
            
        end
        function obj = computeForceVector(obj)
            % Compute force vector of the system
            
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            
            % Obtain force vector. Loop through all bodies and extract the
            % total force for the body.
            forceVector = zeros(3*nBodies,1);
            iF = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    forceVector((3*iF-2:3*iF),1) = obj.myBodies{iB}.myTotalForce;
                    iF = iF + 1;
                end
            end
            obj.myForceVector = forceVector;
        end
        function  obj = computeJpMatrixTotal(obj)
            % Compute the total J matrix for all bodies.
            
            % Determine number of bodies in system and see the J matrix
            if (obj.myBodyIsGround == 1)
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
            end
            JpMatrixTotal = zeros(4*nBodies,4*nBodies);
            
            % Loop through all bodies. Extract Jp for the body. Put Jp into
            % the total matrix.
            iC = 1;
            for iB = 1:obj.myNumBodies
                if (obj.myBodies{iB}.myIsGround == 0)
                    obj.myBodies{iB}.computeJpMatrix();
                    JpMatrix = obj.myBodies{iB}.myJpMatrix;
                    JpMatrixTotal((4*iC-3):4*iC,(4*iC-3):4*iC) = JpMatrix;
                    iC = iC + 1;
                end
            end
            obj.myJpMatrixTotal = JpMatrixTotal;
        end
        function obj = storeSystemState(obj, stepNumber)
            % Store the current state of each body
            % Function inputs:
            %   stepNumber : int
            % Current step number of the analysis
            
            % Extract vectors that contain the state info for all bodies.
            rVec = obj.myR;
            rDotVec = obj.myRDot;
            rDDotVec = obj.myRDDot;
            pVec = obj.myP;
            pDotVec = obj.myPDot;
            pDDotVec = obj.myPDDot;
            time = obj.myTime;
            
            % Add those vectors into the system level total simulation
            % matrix
            obj.myRTotal(:,stepNumber) = rVec;
            obj.myRDotTotal(:,stepNumber) = rDotVec;
            obj.myRDDotTotal(:,stepNumber) = rDDotVec;
            obj.myPTotal(:,stepNumber) = pVec;
            obj.myPDotTotal(:,stepNumber) = pDotVec;
            obj.myPDDotTotal(:,stepNumber) = pDDotVec;
            
            % Reshape the vectors into matrices.
            nBodies = obj.myNumBodies;
            r = reshape(rVec,[3,nBodies]);
            rDot = reshape(rDotVec,[3,nBodies]);
            rDDot = reshape(rDDotVec,[3,nBodies]);
            p = reshape(pVec,[4,nBodies]);
            pDot = reshape(pDotVec,[4,nBodies]);
            pDDot = reshape(pDDotVec,[4,nBodies]);
            
            % Loop through each body. Store the info for each body within
            % its data structure. The resulting data structures will have
            % time and state info across all time steps with each column
            % representing a different time step
            nBodies = obj.myNumBodies;
            for iB = 1:nBodies
                % Store just this time step for the body.
                obj.myBodies{iB}.myTime = time;
                obj.myBodies{iB}.myR = r(:,iB);
                obj.myBodies{iB}.myRDot = rDot(:,iB);
                obj.myBodies{iB}.myRDDot = rDDot(:,iB);
                obj.myBodies{iB}.myP = p(:,iB);
                obj.myBodies{iB}.myPDot = pDot(:,iB);
                obj.myBodies{iB}.myPDDot = pDDot(:,iB);
                
                % Add this to the total time steps for all bodies.
                obj.myBodies{iB}.myTimeTotal(stepNumber) = time;
                obj.myBodies{iB}.myRTotal(:,stepNumber) = r(:,iB);
                obj.myBodies{iB}.myRDotTotal(:,stepNumber) = rDot(:,iB);
                obj.myBodies{iB}.myRDDotTotal(:,stepNumber) = rDDot(:,iB);
                obj.myBodies{iB}.myPTotal(:,stepNumber) = p(:,iB);
                obj.myBodies{iB}.myPDotTotal(:,stepNumber) = pDot(:,iB);
                obj.myBodies{iB}.myPDDotTotal(:,stepNumber) = pDDot(:,iB);
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
            
            % Update the position.
            obj.updateSystemState(rNew, [], [], pNew, [], [], obj.myTime);
            
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
                rInitVec = rInit(4:end,1);
                pInitVec = pInit(5:end,1);
                nBodies = obj.myNumBodies - 1;
            else
                nBodies = obj.myNumBodies;
                rInitVec = rInit;
                pInitVec = pInit;
            end
            
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
                
                % Update the position. If this is the last time through the
                % loop, this will allow us to update the system to the
                % final solution.
                obj.updateSystemState(rNew, [], [], pNew, [], [], time);
                
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
            
            % Update the velocities
            obj.updateSystemState([], rDot, [], [], pDot, [], time);
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
            
            % Update the accelerations
            obj.updateSystemState([], [], rDDot, [], [], pDDot, time);
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
                    obj.computeConstraintProperties(iC, time, phiFlag, 0, 0, 0, 0, 0, 0, 0, 0);
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
                    obj.computeConstraintProperties(iC, time, phiFlag, 0, 0, 0, 0, 0, 0, 0, 0);
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
                    phiP(iP,:) = 0.5*(p'*p) - 0.5;
                    iP = iP + 1;
                end
            end
            
            % Update myPhiP
            obj.myPhiP = phiP;
        end
        function obj = computePhiFullJacobian(obj)
            % Compute the jacobian of phiFull (i.e. partial derivative of
            % phiFull w.r.t. the generalized coordinates)
            
            % Seed the Jacobian. There will be 7 less columns and one less row in the
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
                obj.computeConstraintProperties(iC, time, 0, 0, 0, phiPartialRFlag, phiPartialPFlag, 0, 0, 0, 0);
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
                obj.computeConstraintProperties(iC, time, 0, nuFlag, 0, 0, 0, 0, 0, 0, 0);
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
            
            % Loop through each constraint and compute gamma
            for iC = 1:nKDconst
                gammaFlag = 1;
                obj.computeConstraintProperties(iC, time, 0, 0, gammaFlag, 0, 0, 0, 0, 0, 0);
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
        function obj = updateSystemState(obj, rVec, rDotVec, rDDotVec, pVec, pDotVec, pDDotVec, time)
            % Update position, orientation, and time derivatives of each
            % for each body at a specific time step.
            %
            % Function inputs:
            % rVec : 3N x 1 vector, where N is the number of bodies
            %   vector containing the current 3D positon of each body in
            %   the global reference frame
            %
            % rDotVec : 3N x 1 vector, where N is the number of bodies
            %   vector containing the current time derivative of the 3D
            %   positon of each body in the global reference frame
            %
            % rDDotVec : 3N x 1 vector, where N is the number of bodies
            %   vector containing the current 2nd time derivative of the 3D
            %   positon of each body in the global reference frame
            %
            % pMatrix : 4N x 1 vector, where N is the number of bodies
            %   vector containing the current Euler parameters of each body
            %
            % pDotVec : 4N x 1 vector, where N is the number of bodies
            %   vector containing the current time derivative of the
            %   Euler parameters of each body
            %
            % pDDotVec : 4N x 1 vector, where N is the number of bodies
            %   vector containing the current time derivative of the
            %   Euler parameters of each body
            
            % Update system variables if they are not empty. If they are
            % empty it means those variables are not being updated at this
            % time.
            if ~isempty(rVec)
                obj.myR = rVec;
            end
            if ~isempty(rDotVec)
                obj.myRDot = rDotVec;
            end
            if ~isempty(rDDotVec)
                obj.myRDDot = rDDotVec;
            end
            if ~isempty(pVec)
                obj.myP = pVec;
            end
            if ~isempty(pDotVec)
                obj.myPDot = pDotVec;
            end
            if ~isempty(pDDotVec)
                obj.myPDDot = pDDotVec;
            end
            obj.myTime = time;
            
            % Update individual bodies. If the desired vector is empty,
            % this means that specific value is not being updated at the
            % moment. Just keep this values the same.
            nBodies = obj.myNumBodies;
            for iB = 1:nBodies
                if ~isempty(pVec)
                    pMatrix = reshape(pVec,[4,nBodies]);
                    p = pMatrix(:,iB);
                else
                    p = obj.myBodies{iB}.myP;
                end
                if ~isempty(pDotVec)
                    pDotMatrix = reshape(pDotVec,[4,nBodies]);
                    pDot = pDotMatrix(:,iB);
                else
                    pDot = obj.myBodies{iB}.myPDot;
                end
                if ~isempty(pDDotVec)
                    pDDotMatrix = reshape(pDDotVec,[4,nBodies]);
                    pDDot = pDDotMatrix(:,iB);
                else
                    pDDot = obj.myBodies{iB}.myPDDot;
                end
                if ~isempty(rVec)
                    rMatrix = reshape(rVec,[3,nBodies]);
                    r = rMatrix(:,iB);
                else
                    r = obj.myBodies{iB}.myR;
                end
                if ~isempty(rDotVec)
                    rDotMatrix = reshape(rDotVec,[3,nBodies]);
                    rDot = rDotMatrix(:,iB);
                else
                    rDot = obj.myBodies{iB}.myRDot;
                end
                if ~isempty(rDDotVec)
                    rDDotMatrix = reshape(rDDotVec,[3,nBodies]);
                    rDDot = rDDotMatrix(:,iB);
                else
                    rDDot = obj.myBodies{iB}.myRDDot;
                end
                obj.myBodies{iB}.updateBody(p, pDot, pDDot, r, rDot, rDDot, time);
            end
        end
        function obj = addJoint(obj,jointType,attributes)
            % Add a specific type of joint to the system. Current options
            % are spherical, cylindrical, translational, revolute, and
            % universal.
            
            % All of these joints are kinematic constraints. None are
            % driving constraints.
            %
            % Function inputs:
            % jointType : int
            %   Type of joint. Possible options 'spherical', 'cylindrical',
            %   'translational', 'revolute', or 'universal'.
            %
            % attributes : struct
            %   Structure containing all of the necessary attributes for
            %   your desired joint. This function will check to make
            %   sure all necessary attributes are provided.
            switch jointType
                case 'spherical'
                    obj.createSphericalJoint(attributes);
                case 'cylindrical'
                    obj.createCylindricalJoint(attributes);
                case 'translational'
                    obj.createTranslationalJoint(attributes);
                case 'revolute'
                    obj.createRevoluteJoint(attributes);
                case 'universal'
                    obj.createUniversalJoint(attributes);
                case 'revolute-cylindrical'
                    obj.createRevoluteCylindricalJoint(attributes);
            end
        end
        function obj = createRevoluteCylindricalJoint(obj,attributes)
            % Create a revolute-cylindrical composite joint in this system.
            %
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   translational joint.
            %   body1 = first body in joint. This is the body directly
            %   connected to the revolute joint.
            %   body2 = second body in joint. This is the body directly
            %   connected to the cylindrical joint.
            %   pointOnBody1 = 3D location of the point on body1 used to define the
            %   vector for the b2 constraint in this joint.
            %   pointOnBody2 = 3D location of the point on body2 used to define the
            %   vector for the b2 constraint in this joint.
            %   vectorOnBody1 = Vector on body1 used to define the
            %   revolute joint. This vector is orthogonal to vectorOnBody2.
            %   vector1OnBody2 = Vector on body2 used to define the
            %   translation axis of this body.. This vector is orthogonal to vectorOnBody1.
            %   vector2OnBody2 = 1st vector on body2 used to define the
            %   plane orthogonal to the translational axis. This is
            %   orthogonal to vector1OnBody2.
            %   vector3OnBody2 = 2nd vector on body2 used to define the
            %   plane orthogonal to the translational axis. This is
            %   orthogonal to vector1OnBody2.
            
            necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} ...
                {'pointOnBody2'} {'vectorOnBody1'} {'vector1OnBody2'}  {'vector2OnBody2'}  {'vector3OnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal to zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the DP1 constraint in this joint
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            aBarI = a.vectorOnBody1(:);
            aBarJ = a.vector1OnBody2(:);
            newConstraint1 = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the B2 constraint in this joint
            sBarIP = a.pointOnBody1(:);
            sBarJQ = a.pointOnBody2(:);
            bBarJ = a.vector2OnBody2(:);
            cBarJ = a.vector3OnBody2(:);
            
            % Check to make sure the vectors on body2 are orthogonal
            if (dot(aBarJ,bBarJ)~= 0) || (dot(aBarJ,cBarJ) ~= 0) || (dot(bBarJ,cBarJ)~= 0)
                error('Vectors on body2 must be orthogonal');
            end
            
            % Define DP2 constraints for t
            newConstraint2 = DP2constraint(bodyJ, bodyI, bBarJ, sBarJQ, sBarIP, ft, ftDot, ftDDot, a.constraintName);
            newConstraint3 = DP2constraint(bodyJ, bodyI, cBarJ, sBarJQ, sBarIP, ft, ftDot, ftDDot, a.constraintName);
%             newConstraint2 = DP2constraint(bodyI, bodyJ, bBarJ, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
%             newConstraint3 = DP2constraint(bodyI, bodyJ, cBarJ, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);

            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 3;
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
        end
        function obj = createUniversalJoint(obj,attributes)
            % Create a universal joint in this system.
            %
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   translational joint.
            %   body1 = first body in joint
            %   body2 = second body in joint
            %   pointOnBody1 = 3D location of the point on body1 used to define the
            %   spherical joint portion of the universal joint.
            %   pointOnBody2 = 3D location of the point on body2 used to define the
            %   spherical joint portion of the universal joint.
            %   vectorOnBody1 = Vector on body1 used to define the
            %   universal joint. This vector is orthogonal to vectorOnBody2.
            %   vectorOnBody2 = Vector on body2 used to define the
            %   universal joint. This vector is orthogonal to vectorOnBody1.
            
            necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} ...
                 {'pointOnBody2'} {'vectorOnBody1'} {'vectorOnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the attributes needed for the 3 CD constraints
            % that define the spherical joint portion of the universal
            % joint.
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            sBarIP = a.pointOnBody1;
            sBarJQ = a.pointOnBody2;
            cVec1 = [1 0 0]';
            cVec2 = [0 1 0]';
            cVec3 = [0 0 1]';
            
            % Create the 3 CD constraints for the spherical joint
            newConstraint1 = CDconstraint(bodyI, bodyJ, cVec1, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint2 = CDconstraint(bodyI, bodyJ, cVec2, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint3 = CDconstraint(bodyI, bodyJ, cVec3, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the attributes needed for the DP1 constraint in the
            % universal joint. Create the constraint.
            bodyI = a.body1;
            bodyJ = a.body2;
            aBarI = a.vectorOnBody1;
            aBarJ = a.vectorOnBody2;
            newConstraint4 = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            newConstraint4.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 4;
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
            obj.myConstraints{nConst + 4} = newConstraint4;
        end
        function obj = createTranslationalJoint(obj,attributes)
            % Create a translational joint in this system.
            %
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   translational joint.
            %   body1 = first body in joint
            %   body2 = second body in joint
            %   pointOnBody1 = 3D location of a point along the translational axis of body1
            %   pointOnBody2 = 3D location of a point along the translational axis of body2
            %   vector1OnBody1 = 1st vector to define the plane of the
            %   cylindrical joint on body1
            %   vector2OnBody1 = 2nd vector to define the plane of the
            %   cylindrical joint on body1.
            %   vector1OnBody2 = Vector on body2 that is orthogonal to the
            %   plane on body1.
            %   vector2OnBody2 = Vector on body2 that is parallel to the
            %   plane on body1 and orthogonal to vector1OnBody1.
            
             necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} ...
                 {'pointOnBody2'} {'vector1OnBody1'} {'vector2OnBody1'} {'vector1OnBody2'} {'vector2OnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the attributes needed for the 2 DP1 constraints
            % that define a B1 constraint
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            aBarI = a.vector1OnBody1;
            bBarI = a.vector2OnBody1;
            cBarJ = a.vector1OnBody2;
            newConstraint1 = DP1constraint(bodyI, bodyJ, aBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint2 = DP1constraint(bodyI, bodyJ, bBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the attributes needed for the 2 DP2 constraints
            % that define a B2 constraint
            sBarIP = a.pointOnBody1;
            sBarJQ = a.pointOnBody2;
            newConstraint3 = DP2constraint(bodyI, bodyJ, aBarI, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint4 = DP2constraint(bodyI, bodyJ, bBarI, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the attributes needed for the 3rd DP1 constraint
            aBarJ = a.vector2OnBody2;
            newConstraint5 = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            newConstraint4.myIsKinematic = 1;
            newConstraint5.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 5;
            
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
            obj.myConstraints{nConst + 4} = newConstraint4;
            obj.myConstraints{nConst + 5} = newConstraint5;            
        end
        function obj = createCylindricalJoint(obj,attributes)
            % Create a cylindrical joint in this system.
            %
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   cylindrical joint.
            %   body1 = first body in joint
            %   body2 = second body in joint
            %   pointOnBody1 = 3D location of a point along the translational axis of body1
            %   pointOnBody2 = 3D location of a point along the translational axis of body2
            %   vector1OnBody1 = 1st vector to define the plane of the
            %   cylindrical joint on body1
            %   vector2OnBody1 = 2nd vector to define the plane of the
            %   cylindrical joint on body1.
            %   vectorOnBody2 = Vector on body2 that is orthogonal to the
            %   plane on body1.
            
            necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} {'pointOnBody2'} {'vector1OnBody1'} {'vector2OnBody1'} {'vectorOnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the attributes needed for the 2 DP1 constraints
            % that define a B1 constraint
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            aBarI = a.vector1OnBody1;
            bBarI = a.vector2OnBody1;
            cBarJ = a.vectorOnBody2;
            newConstraint1 = DP1constraint(bodyI, bodyJ, aBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint2 = DP1constraint(bodyI, bodyJ, bBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the attributes needed for the 2 DP2 constraints
            % that define a B2 constraint
            sBarIP = a.pointOnBody1;
            sBarJQ = a.pointOnBody2;
            newConstraint3 = DP2constraint(bodyI, bodyJ, aBarI, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint4 = DP2constraint(bodyI, bodyJ, bBarI, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);

            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            newConstraint4.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 4;

            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
            obj.myConstraints{nConst + 4} = newConstraint4;
        end
        function obj = createRevoluteJoint(obj,attributes)
            % Create a revolute joint in this system.
            %
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   revolute joint.
            %   body1 = first body in joint
            %   body2 = second body in joint
            %   pointOnBody1 = 3D location of the spherical joint on body1
            %   pointOnBody2 = 3D location of the spherical joint on body2
            %   vector1OnBody1 = 1st vector to define the plane of the
            %   revolute joint on body1
            %   vector2OnBody1 = 2nd vector to define the plane of the
            %   revolute joiny on body1.
            %   vectorOnBody2 = Vector on body2 that is orthogonal to the
            %   plane on body1.
            
            necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} {'pointOnBody2'} {'vector1OnBody1'} {'vector2OnBody1'} {'vectorOnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the attributes needed for the 3 CD constraints
            % that define a spherical joint.
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            sBarIP = a.pointOnBody1;
            sBarJQ = a.pointOnBody2;
            cVec1 = [1 0 0]';
            cVec2 = [0 1 0]';
            cVec3 = [0 0 1]';
            
            % Create the 3 CD constraints for the spherical joint
            newConstraint1 = CDconstraint(bodyI, bodyJ, cVec1, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint2 = CDconstraint(bodyI, bodyJ, cVec2, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint3 = CDconstraint(bodyI, bodyJ, cVec3, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            
            % Define the attributes need for the two DP1 constraints
            aBarI = a.vector1OnBody1;
            bBarI = a.vector2OnBody1;
            cBarJ = a.vectorOnBody2;
            newConstraint4 = DP1constraint(bodyI, bodyJ, aBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint5 = DP1constraint(bodyI, bodyJ, bBarI, cBarJ, ft, ftDot, ftDDot, a.constraintName);
            
            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            newConstraint4.myIsKinematic = 1;
            newConstraint5.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 5;
            
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
            obj.myConstraints{nConst + 4} = newConstraint4;
            obj.myConstraints{nConst + 5} = newConstraint5; 
        end
        function obj = createSphericalJoint(obj,attributes)
            % Create a spherical joint in  this system.
            % 
            % Function inputs:
            % attributes : structure
            %   Structure containing the necessary attributes for a
            %   spherical joint.
            %   body1 = first body in joint
            %   body2 = second body in joint
            %   pointOnBody1 = 3D location of the spherical joint on body1
            %   pointOnBody2 = 3D location of the spherical joint on body2
            
            necessaryAttributes = [{'body1'} {'body2'} {'pointOnBody1'} {'pointOnBody2'}];
            
            % Check to make sure attributes are provided.
            for iA = 1:length(necessaryAttributes)
                if ~isfield(attributes,necessaryAttributes{iA})
                    error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                end
            end
            
            % Tell user constraint name is optional if it is not
            % provided
            if ~isfield(attributes,'constraintName')
                disp('constraintName not provided. Setting to default');
                attributes.constraintName = [constraintType ' constraint'];
            end
            
            % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
            % constraint.
            ft = @(t)0;
            ftDot = @(t)0;
            ftDDot = @(t)0;
            
            % Define the attributes needed for the 3 CD constraints
            % that define a spherical joint.
            a = attributes;
            bodyI = a.body1;
            bodyJ = a.body2;
            sBarIP = a.pointOnBody1;
            sBarJQ = a.pointOnBody2;
            cVec1 = [1 0 0]';
            cVec2 = [0 1 0]';
            cVec3 = [0 0 1]';
            
            % Create the 3 CD constraints for the spherical joint
            newConstraint1 = CDconstraint(bodyI, bodyJ, cVec1, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint2 = CDconstraint(bodyI, bodyJ, cVec2, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            newConstraint3 = CDconstraint(bodyI, bodyJ, cVec3, sBarIP, sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            
            
            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = 1;
            newConstraint2.myIsKinematic = 1;
            newConstraint3.myIsKinematic = 1;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if isempty(obj.myNumKinematicConstraints)
                obj.myNumKinematicConstraints = 0;
            end
            obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 3;
            
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with the new constraints
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            obj.myConstraints{nConst + 3} = newConstraint3;
        end
        function obj = addIntermediateConstraint(obj,isKinematic,constraintType,attributes)
            % Add an intermediate constraint to the multibody system. When
            % an intermediate constraint is added it basically creates a
            % number of basic constraints.
            %
            % Function inputs:
            % isKinematic : int
            %   Flag to determine if this is a kinematic or driving
            %   constraint. 1 if constraint is kinematic. 0 if constraint
            %   is driving
            %
            % constraintType : int
            %   Type of intermediate constraint. Possible options 'b1' or
            %   'b2'
            %
            % attributes : struct
            %   Structure containing all of the necessary attributes for
            %   your desired constraint. This function will check to make
            %   sure all necessary attributes are provided.
            
            switch constraintType
                case 'b1'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'aBarI'}  {'bBarI'} {'cBarJ'}];
                    
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
              
                    % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
                    % constraint.
                    ft = 0;
                    ftDot = 0;
                    ftDDot = 0;
                    
                    % Add the first DP1 constraint between two bodies in the system
                    a = attributes;
                    newConstraint1 = DP1constraint(a.bodyI, a.bodyJ, a.aBarI, a.cBarJ, ft, ftDot, ftDDot, a.constraintName);
                    newConstraint2 = DP1constraint(a.bodyI, a.bodyJ, a.bBarI, a.cBarJ, ft, ftDot, ftDDot, a.constraintName);
                    
                case 'b2'
                    % Check to make sure all necessary attributes have been
                    % provided
                    necessaryAttributes = [{'bodyI'} {'bodyJ'} {'aBarI'}  {'bBarI'} {'sBarIP'} {'sBarJQ'}];
                    
                    for iA = 1:length(necessaryAttributes)
                        if ~isfield(attributes,necessaryAttributes{iA})
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        disp('constraintName not provided. Setting to default');
                        attributes.constraintName = [constraintType ' constraint'];
                    end
              
                    % Set f(t) fDot(t) and fDDot(t) equal tro zero for this
                    % constraint.
                    ft = 0;
                    ftDot = 0;
                    ftDDot = 0;
                    
                    % Add the first DP1 constraint between two bodies in the system
                    a = attributes;
                    newConstraint1 = DP2constraint(a.bodyI, a.bodyJ, a.aBarI, a.sBarIP, a.sBarJQ, ft, ftDot, ftDDot, a.constraintName);
                    newConstraint2 = DP2constraint(a.bodyI, a.bodyJ, a.bBarI, a.sBarIP, a.sBarJQ, ft, ftDot, ftDDot, a.constraintName);
            end
            
            % Set flag for kinematic vs driving constraint
            newConstraint1.myIsKinematic = isKinematic;
            newConstraint2.myIsKinematic = isKinematic;
            
            % Update count of kinematic constraints and number of driving
            % constraints
            if (isKinematic)
                if isempty(obj.myNumKinematicConstraints)
                    obj.myNumKinematicConstraints = 0;
                end
                obj.myNumKinematicConstraints = obj.myNumKinematicConstraints + 2;
            else
                if isempty(obj.myNumDrivingConstraints)
                    obj.myNumDrivingConstraints = 0;
                end
                obj.myNumDrivingConstraints = obj.myNumDrivingConstraints + 2;
            end
            
            % Current number of constraints
            nConst = obj.myNumConstraints;
            
            % Update system with both new constraint
            obj.myConstraints{nConst + 1} = newConstraint1;
            obj.myConstraints{nConst + 2} = newConstraint2;
            
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
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
                        end
                    end
                    
                    % Tell user constraint name is optional if it is not
                    % provided
                    if ~isfield(attributes,'constraintName')
                        error('constraintName not provided. Setting to default');
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
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
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
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
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
                            error(['ERROR: Must provide ' necessaryAttributes{iA} ' for ' constraintType ' constraint.']);
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
        function obj = computeConstraintProperties(obj, constraintNumber, time, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag, constraintForcePartialRFlag, constraintForcePartialPFlag, constraintTorquePartialRFlag, constraintTorquePartialPFlag)
            % Determine type of constraint
            cType = obj.myConstraints{constraintNumber}.myConstraintType;
            
            % Call the correct constraint class based off constraint type.
            switch cType
                case 'DP1'
                    obj.myConstraints{constraintNumber}.computeDP1constraint(obj, time, phiFlag, nuFlag, gammaFlag, ...
                        phiPartialRFlag, phiPartialPFlag, constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                        constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber);
                    
                case 'DP2'
                    obj.myConstraints{constraintNumber}.computeDP2constraint(obj, time, phiFlag, nuFlag, gammaFlag, ...
                        phiPartialRFlag, phiPartialPFlag, constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                        constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber);
                    
                case 'D'
                    obj.myConstraints{constraintNumber}.computeDconstraint(obj, time, phiFlag, nuFlag, gammaFlag, ...
                        phiPartialRFlag, phiPartialPFlag, constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                        constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber);
                    
                case 'CD'
                    obj.myConstraints{constraintNumber}.computeCDconstraint(obj, time, phiFlag, nuFlag, gammaFlag, ...
                        phiPartialRFlag, phiPartialPFlag, constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                        constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber);
            end
        end
        
    end
    % Methods block with no attributes
    methods
        function myNumBodies = get.myNumBodies(obj)
            % Calculate number of bodies in system
            myNumBodies = length(obj.myBodies);
        end
        function myNumBodiesMinusGround = get.myNumBodiesMinusGround(obj)
            % Calculate number of bodies in system minus the ground. If
            % there is no ground in the system, myNumBodiesMinusGround =
            % myNumBodies
            if (obj.myBodyIsGround == 1)
                myNumBodiesMinusGround = length(obj.myBodies) - 1;
            else
                myNumBodiesMinusGround = length(obj.myBodies);
            end
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
        function myNumTimeSteps = get.myNumTimeSteps(obj)
            % Compute current number of time steps completed for simulation
            myNumTimeSteps = length(obj.myTimeTotal);
        end
    end
    
end