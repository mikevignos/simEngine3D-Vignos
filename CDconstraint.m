classdef CDconstraint < handle
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'CD'; % String. Type of basic constraint. Currently only DP1 or CD
        myBodyI; % First body in constraint
        myBodyJ; % Second body in constraint
        myCoordVec; % Coordinate that is being constraint (e.g. for x constraint, myCoordVec = [1 0 0]')
        mysBarIP; % Point on body I being constrained in body I reference frame
        mysBarJQ; % Point on body J being constrained in body J reference frame
        myFt; % Function representing value of constraint
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
        function obj = CDconstraint(bodyI, bodyJ, coordVec, sBarIP, sBarJQ, ft, ftDot, ftDDot, constraintName)
            % Store name sent in by user
            if nargin < 9
                obj.myConstraintName = 'CD Constraint';
            else
                obj.myConstraintName = constraintName;
            end
            
            % Store attributes.
            obj.myBodyI = bodyI;
            obj.myBodyJ = bodyJ;
            obj.myCoordVec = coordVec;
            obj.mysBarIP = sBarIP;
            obj.mysBarJQ = sBarJQ;
            obj.myFt = ft;
            obj.myFtDot = ftDot;
            obj.myFtDDot = ftDDot;
        end
        function obj = computeCDconstraint(obj, sys, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
            % Computes necessary quantities for CD constraint
            % Possible quantities to compute are described in the
            % properties of this class.
            %
            % Function inputs:
            % allBodies : struct
            %   Data structure containing all of the bodies in the system.
            %   Each body in this structure should be defined using
            %   the body class.
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
        end
        function obj = computePhi(obj,sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            ft = obj.myFt;
            t = obj.myTime;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % Compute necessary parameters for body I. If bodyI is
            % the ground. Hardcode in parameters for this case.
            % Orientation matrix
            if (isGroundI == 1)
                Ai = eye(3,3);
                ri = zeros(3,1);
            else
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                ri = sys.myBodies{bodyI}.myR;
            end
            
            % Compute necessary parameters for body J. If bodyJ is
            % the ground. Hardcode in parameters for this case.
            if (isGroundJ == 1)
                Aj = eye(3,3);
                rj = zeros(3,1);
            else
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
                rj = sys.myBodies{bodyJ}.myR;
            end
            
            % Compute phi
            ftVal = ft(t);
            phi = coordVec'*(rj + Aj*sBarJQ - ri - Ai*sBarIP) - ftVal;
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
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            t = obj.myTime;
            ftDDot = obj.myFtDDot;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % Compute necessary parameters for body I. Check for special
            % case when bodyI is the ground.
            if (isGroundI == 1)
                pDotI = zeros(4,1);
                BdotI = zeros(3,4);
            else
                % Time derivatives of Euler parameters
                pDotI = sys.myBodies{bodyI}.myPDot;
                
                % Time derivative of B matrix
                sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeBDot(sBarIP);
                BdotI = sys.myBodies{bodyI}.myBDot;
            end
            
            % Compute necessary parameters for body J. Check for special
            % case when bodyJ is the ground.
            if (isGroundJ == 1)
                pDotJ = zeros(4,1);
                BdotJ = zeros(3,4);
            else
                % Time derivatives of Euler parameters
                pDotJ = sys.myBodies{bodyJ}.myPDot;
                
                % Compute B(pDot, aBar)
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeBDot(sBarJQ);
                BdotJ = sys.myBodies{bodyJ}.myBDot;
            end
            
            %Compute right hand side of acceleration equation
            ftDDotVal = ftDDot(t);
            gamma = coordVec'*BdotI*pDotI - coordVec'*BdotJ*pDotJ + ftDDotVal;
            obj.myGamma = gamma;
        end
        function obj = computePhiPartialR(obj, sys)
            % Extract needed attributes
            bodyJ = obj.myBodyJ;
            bodyI = obj.myBodyI;
            coordVec = obj.myCoordVec;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1)
                phiPartialRI = -coordVec';
                phiPartialR = phiPartialRI;
                
                % If bodyI is the ground, this body contributes
                % nothing to the Jacobian
            elseif (isGroundI == 1)
                phiPartialRJ = coordVec';
                phiPartialR = phiPartialRJ;
            else
                phiPartialRI = -coordVec';
                phiPartialRJ = coordVec';
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
            coordVec = obj.myCoordVec;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute B(p, sBar) for bodyI
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                % Compute phiPartialP
                phiPartialPI = -coordVec'*BmatrixI;
                phiPartialP = phiPartialPI;
                
                % If bodyI is the ground, this body contributes
            % nothing to the Jacobian
            elseif (isGroundI == 1);
                % Compute B(p, sBar) for bodyJ
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                
                % Compute phiPartialP
                phiPartialPJ = coordVec'*BmatrixJ;
                phiPartialP = phiPartialPJ;
            else
                % Compute B(p, sBar) for bodyI and body J
                sys.myBodies{bodyI}.computeB(sBarIP);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyJ}.computeB(sBarJQ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                
                % Compute phiPartialP
                phiPartialPI = -coordVec'*BmatrixI;
                phiPartialPJ = coordVec'*BmatrixJ;
                phiPartialP = [phiPartialPI phiPartialPJ];
            end
            
            obj.myPhiPartialP = phiPartialP;
        end
    end
    
end
