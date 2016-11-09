classdef Dconstraint < handle
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'D'; % String. Type of basic constraint. Currently only DP1 or CD
        myIsKinematic; % Flag for if the constraint is a kinematic (1) or driving (0) constraint
        myBodyI; % First body in constraint
        myBodyJ; % Second body in constraint
        mysBarIP; % Point on body I being constrained in body I reference frame
        mysBarJQ; % Point on body J being constrained in body J reference frame
        myFt; % Function representing the square of the distance between points mysBarIP and mysBarJQ
        myFtDot; % Function representing first time derivative of f(t)
        myFtDDot; % Function representing second time derivative of f(t)
        myTime; % Value of current time step
        
        % Properties of constraint that can be computed
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates (i.e. r)
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates (i.e. p)
    end
    
    methods
        function obj = Dconstraint(bodyI, bodyJ, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName)
            % Store name sent in by user
            if nargin < 8
                obj.myConstraintName = 'D Constraint';
            else
                obj.myConstraintName = constraintName;
            end
            
            % Store attributes.
            obj.myBodyI = bodyI;
            obj.myBodyJ = bodyJ;
            obj.mysBarIP = sBarIP;
            obj.mysBarJQ = sBarJQ;
            obj.myFt = ft;
            obj.myFtDot = ftDot;
            obj.myFtDDot = ftDDot;
        end
        function obj = computeDconstraint(obj, sys, t, phiFlag, nuFlag, gammaFlag, ...
                phiPartialRFlag, phiPartialPFlag, ...
                constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber)
            % Computes necessary quantities for D constraint
            % Possible quantities to compute are described in the
            % properties of this class.
            %
            % Function inputs:
            % sys : class
            %   multibodySystem class. This class contains all of the info
            %   about the bodies in your multibody system.
            %
            % t : double
            %   Value of current time step
            %
            % phiFlag : int
            %   Flag indicating if phi should be computed
            %
            % nuFlag : int
            %   Flag indicating if nu should be computed
            %
            % gammaFlag : int
            %   Flag indicating if gamma should be computed
            %
            % phiPartialRFlag : int
            %   Flag indicating if phiPartialR should be computed
            %
            % phiPartialPFlag : int
            %     Flag indicating if phiPartialP should be computed
            
            % Store time
            obj.myTime = t;
            
            % Compute desired parameters
            if (phiFlag == 1)
                obj = obj.computePhi(sys);
            end
            if (nuFlag == 1)
                obj = obj.computeNu();
            end
            if (gammaFlag == 1)
                obj = obj.computeGamma(sys);
            end
            if (phiPartialRFlag == 1)
                obj = obj.computePhiPartialR(sys);
            end
            if (phiPartialPFlag == 1)
                obj = obj.computePhiPartialP(sys);
            end
            if (constraintForcePartialRFlag == 1)
                obj.computeConstraintForcePartialR(sys, constraintNumber);
            end
            if (constraintForcePartialPFlag == 1)
                obj.computeConstraintForcePartialP(sys, constraintNumber);
            end
            if (constraintTorquePartialRFlag == 1)
                obj.computeConstraintTorquePartialR(sys, constraintNumber);
            end
            if (constraintTorquePartialPFlag == 1)
                obj.computeConstraintTorquePartialP(sys, constraintNumber);
            end
        end
        function obj = computePhi(obj,sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            ft = obj.myFt;
            t = obj.myTime;

            % Compute dij
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % f(t) already represents the square of the distance as a 
            % function of time
            ftVal = ft(t);
            
            % Compute phi
            phi = dij'*dij - ftVal;
            obj.myPhi = phi;
        end
        function obj = computeNu(obj)
            % Extract needed attributes
            ftDot = obj.myFtDot;
            t = obj.myTime;
            
            % Compute nu.
            nu = ftDot(t);
            obj.myNu = nu;
        end
        function obj = computeGamma(obj, sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            t = obj.myTime;
            ftDDot = obj.myFtDDot;
            
            % Compute necessary parameters
            % Time derivatives of Euler parameters
            pDotI = sys.myBodies{bodyI}.myPDot;
            pDotJ = sys.myBodies{bodyJ}.myPDot;
            
            % Time derivative of B matrix
            sys.myBodies{bodyI}.computeBDot(sBarIP);
            BdotSBarIP = sys.myBodies{bodyI}.myBDot;
            
            sys.myBodies{bodyJ}.computeBDot(sBarJQ);
            BdotSBarJQ = sys.myBodies{bodyJ}.myBDot;
            
            % Vector and time derivative of vector
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            dijDot = simEngine3DUtilities.computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ);

            %Compute right hand side of acceleration equation
            ftDDotVal = ftDDot(t);
            gamma = -2*dij'*BdotSBarJQ*pDotJ + 2*dij'*BdotSBarIP*pDotI - 2*(dijDot'*dijDot) + ftDDotVal;
            obj.myGamma = gamma;
        end
        function obj = computePhiPartialR(obj, sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % Compute dij because you need it for all cases
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1)
                phiPartialRI = -2*dij';
                phiPartialR = phiPartialRI;
                
                % If bodyI is the ground, this body contributes
                % nothing to the Jacobian
            elseif (isGroundI == 1)
                phiPartialRJ = 2*dij';
                phiPartialR = phiPartialRJ;
            else
                phiPartialRI = -2*dij';
                phiPartialRJ = 2*dij';
                phiPartialR = [phiPartialRI phiPartialRJ];
            end
            obj.myPhiPartialR = phiPartialR;
        end
        function obj = computePhiPartialP(obj,sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % Compute dij because you need it for all cases
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute B(p, sBar) for bodyI
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSBarIP = sys.myBodies{bodyI}.myB;
                
                % Compute phiPartialP
                phiPartialPI = -2*dij'*BmatrixSBarIP;
                phiPartialP = phiPartialPI;
                
                % If bodyI is the ground, this body contributes
                % nothing to the Jacobian
            elseif (isGroundI == 1);
                % Compute B(p, sBar) for bodyJ
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixSBarJQ = sys.myBodies{bodyJ}.myB;
                
                % Compute phiPartialP
                phiPartialPJ = 2*dij'*BmatrixSBarJQ;
                phiPartialP = phiPartialPJ;
            else
                % Compute B(p, sBar) for bodyI and body J
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSBarIP = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixSBarJQ = sys.myBodies{bodyJ}.myB;
                
                % Compute phiPartialP
                phiPartialPI = -2*dij'*BmatrixSBarIP;
                phiPartialPJ = 2*dij'*BmatrixSBarJQ;
                phiPartialP = [phiPartialPI phiPartialPJ];
            end
            
            obj.myPhiPartialP = phiPartialP;
        end
        function obj = computeConstraintForcePartialR(obj, sys, constraintNumber)
            % Compute the partial derivative of the forces due to this
            % constraint, w.r.t position (r)
            
            % Extract the lagrange multiplier that relates to this
            % constraint
            lambda = sys.myConstraintLagrangeMultipliers(constraintNumber);
            
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                constForcePartialRI = [1 0 0; 0 1 0; 0 0 1];
                
                constForcePartialR = 2*lambda*constForcePartialRI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                constForcePartialRJ = [1 0 0; 0 1 0; 0 0 1];
                
                constForcePartialR = 2*lambda*constForcePartialRJ;
            else
                constForcePartialRI = zeros(6,3);
                constForcePartialRI(1:3,:) = [1 0 0; 0 1 0; 0 0 1];
                constForcePartialRI(4:6,:) = -[1 0 0; 0 1 0; 0 0 1];
                
                constForcePartialRJ = zeros(6,3);
                constForcePartialRJ(1:3,:) = -[1 0 0; 0 1 0; 0 0 1];
                constForcePartialRJ(4:6,:) = [1 0 0; 0 1 0; 0 0 1];
                
                constForcePartialR = 2*lambda*[constForcePartialRI constForcePartialRJ];                
            end
            
            obj.myConstraintForcePartialR = constForcePartialR;
        end
        function obj = computeConstraintForcePartialP(obj, sys, constraintNumber)
            % Compute the partial derivative of the forces due to this
            % constraint, w.r.t orientation (p)
            
            % Extract the lagrange multiplier that relates to this
            % constraint
            lambda = sys.myConstraintLagrangeMultipliers(constraintNumber);
            
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute B(p, sBarIP) for bodyI.
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                constForcePartialPI = BmatrixI;
                
                constForcePartialP = 2*lambda*constForcePartialPI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                % Compute B(p, sBarJQ) for bodyJ.
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                constForcePartialPJ = BmatrixJ;
                
                constForcePartialP = 2*lambda*constForcePartialPJ;
            else
                % Compute B(p, sBarIP) for bodyI.
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                % Compute B(p, sBarJQ) for bodyJ.
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                
                constForcePartialPI = zeros(6,4);
                constForcePartialPI(1:3,:) = BmatrixI;
                constForcePartialPI(4:6,:) = -BmatrixI;
                
                constForcePartialPJ = zeros(6,4);
                constForcePartialPJ(1:3,:) = -BmatrixJ;
                constForcePartialPJ(4:6,:) = BmatrixJ;

                constForcePartialP = 2*lambda*[constForcePartialPI constForcePartialPJ];                
            end
            
            obj.myConstraintForcePartialP = constForcePartialP;
        end
        function obj = computeConstraintTorquePartialR(obj, sys, constraintNumber)
            % Compute the partial derivative of the torques due to this
            % constraint, w.r.t position (r)
            
            % Extract the lagrange multiplier that relates to this
            % constraint
            lambda = sys.myConstraintLagrangeMultipliers(constraintNumber);
            
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
               % Compute B(p, sBarIP) for bodyI.
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                constTorquePartialRI = BmatrixI';
                               
                constTorquePartialR = 2*lambda*constTorquePartialRI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                % Compute B(p, sBarJQ) for bodyJ.
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyI}.myB;
                constTorquePartialRJ = BmatrixJ;
                
                constTorquePartialR = 2*lambda*constTorquePartialRJ;
            else
                % Compute B(p, sBarIP) for bodyI.
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                % Compute B(p, sBarJQ) for bodyJ.
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyI}.myB;
                
                constTorquePartialRI = zeros(8,3);
                constTorquePartialRI(1:4,:) = BmatrixI';
                constTorquePartialRI(5:8,:) = -BmatrixJ;
                
                constTorquePartialRJ = zeros(8,3);
                constTorquePartialRJ(1:4,:) = -BmatrixI';
                constTorquePartialRJ(5:8,:) = BmatrixJ;
                
                constTorquePartialR = 2*lambda*[constTorquePartialRI constTorquePartialRJ];                
            end
            
            obj.myConstraintTorquePartialR = constTorquePartialR;          
        end
        %%%%%%%%%%%STOPPED HERE!!!!!%%%%%%%%%%%%%%%%
        function     obj = computeConstraintTorquePartialP(obj, sys, constraintNumber)
            % Compute the partial derivative of the torques due to this
            % constraint, w.r.t orientation (p)
            
            % Extract the lagrange multiplier that relates to this
            % constraint
            lambda = sys.myConstraintLagrangeMultipliers(constraintNumber);
            
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            aBarI = obj.myaBarI;
            
            % Compute distance vector between the points
            dij = computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute distance vector between the two points, dij
                dij = simEngine3DUtilities.computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                
                % Compute orientation matrix, A, for bodyI
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                               
                % Compute a for bodyJ
                ai = Ai*aBarI;
                
                % Compute the two K matrices that are needed.
                K_aBarI_dij = simEngine3DUtilities.computeKMatrix(aBarI, dij);
                K_sBarIP_ai = simEngine3DUtilities.computeKMatrix(sBarIP, ai);
                
                % Compute B matrices that are needed.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixAbarI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSbarIP = sys.myBodies{bodyI}.myB;
                
                % Compute the X term
                constTorquePartialPI = K_aBarI_dij - K_sBarIP_ai - BmatrixAbarI'*BmatrixSbarIP - BmatrixSbarIP'*BmatrixAbarI;

                constTorquePartialP = 2*lambda*constTorquePartialPI;
                
                % Special case when body I is ground
            elseif (isGroundI == 1);
                % Compute orientation matrix, A, for bodyI
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                               
                % Compute a for bodyJ
                ai = Ai*aBarI;
                
                % Compute the K matrix that is needed.
                K_sBarJQ_ai = simEngine3DUtilities.computeKMatrix(sBarJQ, ai);
                
                % Compute constTorquePartialP
                constTorquePartialPJ = K_sBarJQ_ai;
                constTorquePartialP = 2*lambda*constTorquePartialPJ;
                
            else
                % Compute distance vector between the two points, dij
                dij = simEngine3DUtilities.computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                
                % Compute orientation matrix, A, for bodyI
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                               
                % Compute a for bodyJ
                ai = Ai*aBarI;
                
                % Compute three two K matrices that are needed.
                K_aBarI_dij = simEngine3DUtilities.computeKMatrix(aBarI, dij);
                K_sBarIP_ai = simEngine3DUtilities.computeKMatrix(sBarIP, ai);
                K_sBarJQ_ai = simEngine3DUtilities.computeKMatrix(sBarJQ, ai);
                
                % Compute B matrices that are needed.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixAbarI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSbarIP = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixSbarJQ = sys.myBodies{bodyI}.myB;
                
                % Compute the X term
                Xterm = K_aBarI_dij - K_sBarIP_ai - BmatrixAbarI'*BmatrixSbarIP - BmatrixSbarIP'*BmatrixAbarI;

                % Compute constTorquePartialP
                constTorquePartialPI = zeros(8,4);
                constTorquePartialPI(1:4,:) = Xterm;
                constTorquePartialPI(5:8,:) = BmatrixSbarJQ'*BmatrixAbarI;
                
                constTorquePartialPJ = zeros(8,4);
                constTorquePartialPJ(1:4,:) = BmatrixAbarI'*BmatrixSbarJQ;
                constTorquePartialPJ(5:8,:) = K_sBarJQ_ai;
                
                constTorquePartialP = 2*lambda*[constTorquePartialPI constTorquePartialPJ];
            end
            
            obj.myConstraintTorquePartialP = constTorquePartialP;
        end
    end
    
end
