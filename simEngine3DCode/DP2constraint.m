classdef DP2constraint < handle
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP2'; % String. Type of basic constraint. Currently only DP1 or CD
        myIsKinematic; % Flag for if the constraint is a kinematic (1) or driving (0) constraint
        myBodyI; % First body in constraint
        myBodyJ; % Second body in constraint
        myaBarI; % Vector defined on body I in body I reference frame
        mysBarIP; % Point on body I in body I reference frame
        mysBarJQ; % Point on body J in body J reference frame
        myFt; % Function representing the value of the dot product between myaBarI and the vector between mysBarIP and mysBarJQ.
        myFtDot; % Function representing first time derivative of f(t)
        myFtDDot; % Function representing second time derivative of f(t)
        myTime; % Value of current time step
        
        % Properties of constraint that can be computed
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates (i.e. r)
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates (i.e. p)
        myConstraintForcePartialR; % Partial derivative of this constraint force w.r.t R.
        myConstraintForcePartialP; % Partial derivative of this constraint force w.r.t P.
        myConstraintTorquePartialR; % Partial derivative of this constraint torque w.r.t R.
        myConstraintTorquePartialP; % Partial derivative of this constraint torque w.r.t P.
    end
    
    methods
        function obj = DP2constraint(bodyI, bodyJ, aBarI, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName)
            % Store name sent in by user
            if nargin < 9
                obj.myConstraintName = 'DP2 Constraint';
            else
                obj.myConstraintName = constraintName;
            end
            
            % Make sure all vectors are column vectors
            aBarI  = aBarI(:);
            sBarIP = sBarIP(:);
            sBarJQ = sBarJQ(:);
            
            % Store attributes.
            obj.myBodyI = bodyI;
            obj.myBodyJ = bodyJ;
            obj.myaBarI = aBarI;
            obj.mysBarIP = sBarIP;
            obj.mysBarJQ = sBarJQ;
            obj.myFt = ft;
            obj.myFtDot = ftDot;
            obj.myFtDDot = ftDDot;
            
        end
        function obj = computeDP2constraint(obj, sys, t, phiFlag, nuFlag, gammaFlag, ...
                phiPartialRFlag, phiPartialPFlag, ...
                constraintForcePartialRFlag, constraintForcePartialPFlag, ...
                constraintTorquePartialRFlag, constraintTorquePartialPFlag, constraintNumber)
            % Computes necessary quantities for DP1 constraint
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
            if(nuFlag == 1)
                obj = obj.computeNu();
            end
            if(gammaFlag == 1)
                obj = obj.computeGamma(sys);
            end
            if(phiPartialRFlag == 1)
                obj = obj.computePhiPartialR(sys);
            end
            if(phiPartialPFlag == 1)
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
            aBarI = obj.myaBarI;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            ft = obj.myFt;
            t = obj.myTime;
            
            % Compute necessary parameters
            % Orientation matrix
            sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            % Compute dij
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % Compute phi
            ftVal = ft(t);
            phi = aBarI'*Ai'*dij - ftVal;
            obj.myPhi = phi;
        end
        function obj = computeNu(obj)
            % Extract needed attributes
            ftDot = obj.myFtDot;
            t = obj.myTime;
            
            % Compute nu
            nu = ftDot(t);
            obj.myNu = nu;
        end
        function obj = computeGamma(obj, sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            aBarI = obj.myaBarI;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            ftDDot = obj.myFtDDot;
            t = obj.myTime;
            
            % Compute dij and time derivative of dij
            dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            dijDot = simEngine3DUtilities.computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ);
            
            % Compute necessary parameters
            % Orientation matrix
            sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            % Time derivatives of Euler parameters
            pDotI = sys.myBodies{bodyI}.myPDot;
             pDotJ = sys.myBodies{bodyJ}.myPDot;
             
            % B matrix
            sys.myBodies{bodyI}.computeB(aBarI);
            BmatrixABarI = sys.myBodies{bodyI}.myB;
            
            % Time derivative of B matrix
            sys.myBodies{bodyI}.computeBDot(sBarIP);
            BdotSBarIP = sys.myBodies{bodyI}.myBDot;
            
            sys.myBodies{bodyJ}.computeBDot(sBarJQ);
            BdotSBarJQ = sys.myBodies{bodyJ}.myBDot;
            
            sys.myBodies{bodyI}.computeBDot(aBarI);
            BdotABarI= sys.myBodies{bodyI}.myBDot;         
            
            % a and aDot
            ai = Ai*aBarI;
            aDotI = BmatrixABarI*pDotI;

            %Compute right hand side of acceleration equation
            gamma = -ai'*BdotSBarJQ*pDotJ + ai'*BdotSBarIP*pDotI - dij'*BdotABarI*pDotI - 2*aDotI'*dijDot + ftDDot(t);
            obj.myGamma = gamma;
        end
        function obj = computePhiPartialR(obj, sys)
            % Extract necessary parameters
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            aBarI = obj.myaBarI;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % Compute necessary properties
            % Orientation matrix
            sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            ai = Ai*aBarI;
            
            % Compute phiPartialR
            if (isGroundJ == 1)
                phiPartialRI = -ai';
                phiPartialR = phiPartialRI;
            elseif (isGroundI == 1)
                phiPartialRJ = ai';
                phiPartialR = phiPartialRJ;
            else
                phiPartialRI = -ai';
                phiPartialRJ = ai';
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
            aBarI = obj.myaBarI;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute necessary parameters
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                ai = Ai*aBarI;
                dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                
                % B matrices
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixABarI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSBarIP = sys.myBodies{bodyI}.myB;
                
                % Compute phiPartialP
                phiPartialPI = dij'*BmatrixABarI - ai'*BmatrixSBarIP;
                phiPartialP = phiPartialPI;
                
                 % Special case when body I is ground
            elseif (isGroundI == 1);
                % Compute necessary parameters
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                ai = Ai*aBarI;
                
                % B matrix
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixSBarJQ = sys.myBodies{bodyJ}.myB;
                
                % Compute phiPartialP
                phiPartialPJ = ai'*BmatrixSBarJQ;
                phiPartialP = phiPartialPJ;
            else
                % Compute necessary parameters
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                ai = Ai*aBarI;
                dij = simEngine3DUtilities.computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ);
                
                % B matrices
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixABarI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixSBarIP = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixSBarJQ = sys.myBodies{bodyJ}.myB;

                % Compute phiPartialP
                phiPartialPI = dij'*BmatrixABarI - ai'*BmatrixSBarIP;
                phiPartialPJ = ai'*BmatrixSBarJQ;
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
                constForcePartialRI = zeros(3,3);
                
                constForcePartialR = lambda*constForcePartialRI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                constForcePartialRJ = zeros(3,3);
                
                constForcePartialR = lambda*constForcePartialRJ;
            else
                constForcePartialRI = zeros(6,3);
                constForcePartialRJ = zeros(6,3);
                constForcePartialR = lambda*[constForcePartialRI constForcePartialRJ];                
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
            aBarI = obj.myaBarI;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute B(p, aBar) for bodyI.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;
                constForcePartialPI = -BmatrixI;
                
                constForcePartialP = lambda*constForcePartialPI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                constForcePartialPJ = zeros(3,4);
                
                constForcePartialP = lambda*constForcePartialPJ;
            else
                % Compute B(p, aBar) for bodyI.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;                
                
                constForcePartialPI = zeros(6,4);
                constForcePartialPI(1:3,:) = -BmatrixI;
                constForcePartialPI(4:6,:) = BmatrixI;
                
                constForcePartialPJ = zeros(6,4);

                constForcePartialP = lambda*[constForcePartialPI constForcePartialPJ];                
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
               % Compute B(p, aBar) for bodyI.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;
                constTorquePartialRI = -BmatrixI';
                               
                constTorquePartialR = lambda*constTorquePartialRI;
               
                % Special case when body I is ground
            elseif (isGroundI == 1);
                constTorquePartialRJ = zeros(4,3);
                
                constTorquePartialR = lambda*constTorquePartialRJ;
            else
                % Compute B(p, aBar) for bodyI.
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;                
                
                constTorquePartialRI = zeros(8,3);
                constTorquePartialRI(1:4,:) = -BmatrixI';
                
                constTorquePartialRJ = zeros(8,3);
                constTorquePartialRJ(1:3,:) = BmatrixI';
                
                constTorquePartialR = lambda*[constTorquePartialRI constTorquePartialRJ];                
            end
            
            obj.myConstraintTorquePartialR = constTorquePartialR;          
        end
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

                constTorquePartialP = lambda*constTorquePartialPI;
                
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
                constTorquePartialP = lambda*constTorquePartialPJ;
                
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
                
                constTorquePartialP = lambda*[constTorquePartialPI constTorquePartialPJ];
            end
            
            obj.myConstraintTorquePartialP = constTorquePartialP;
        end
    end
    
end
