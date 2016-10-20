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
        myPhiK; % Vector of kinematic constraints
        myPhiD; % Vector of driving constraints
        myPhiP; % Vector of Euler parameter normalization constraints
        myPhiFull; % Full constraint matrix. Made up of combination of kinematic, driving, and Euler parameter normalization constraints.
        myPhiFullJacobian; % Jacobian of full constraint matrix.
        myNu; % RHS of velocity equation for entire multibody system
        myGamma; % RHS of acceleration equation for entire multibody system
        myTime; % Current time within system.
    end
    
    properties (Dependent)
        myNumBodies; % Number of bodies in the system
        myNumConstraints; % Number of kinematic and driving constraints in the system
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
            newBody = body(bodyNumber, bodyType, isGround, mass, bodyLength);
            obj.myBodies{bodyNumber} = newBody;
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
            
            for iC = 1:nConst
                if (obj.myConstraints{iC}.myIsKinematic == 0)
                    time = obj.myTime;
                    phiFlag = 1;
                    obj.computeConstraintProperties(iC, time, phiFlag, 0, 0, 0, 0);
                    phiD(iC,:) = obj.myConstraints{iC}.myPhi;
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
                phiPartialR = obj.myConstraints(iC).myPhiPartialR;
                phiPartialP = obj.myConstraints(iC).myPhiPartialP;
                
                % If one of the bodies in the system is the ground adjust
                % the body numbers accordingly.
                if (obj.myBodyIsGround == 1)
                    if (obj.myBodies{bodyI}.myIsGround == 1)
                        % Only care about bodyJ in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. This only works because I
                        % am forcing the ground to be body 1.
                        phiFullPartialR(iC,(3*(bodyJ-1) - 2):3*(bodyJ-1)) = phiPartialR;
                        phiFullPartialP(iC,(4*(bodyJ-1) - 2):4*(bodyJ-1)) = phiPartialP;
                        
                    elseif (obj.myBodies{bodyJ}.myIsGround == 1)
                        % Only care about bodyI in this case.
                        % Need to decrease body number by 1 to properly
                        % populate the Jacobian. 
                        phiFullPartialR(iC,(3*(bodyI-1) - 2):3*(bodyI-1)) = phiPartialR;
                        phiFullPartialP(iC,(4*(bodyI-1) - 2):4*(bodyI-1)) = phiPartialP;
                        
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
            
            nBodies = obj.myNumBodies;
            for iB = 1:nBodies
                p = pMatrix(:,iB);
                pDot = pDotMatrix(:,iB);
                r = rMatrix(:,iB);
                rDot = rDotMatrix(:,iB);
                obj.myBodies{iB}.updateBody(p, pDot, r, rDot, time);
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
    end
    
end