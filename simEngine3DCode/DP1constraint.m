classdef DP1constraint < handle
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
        myIsKinematic; % Flag for if the constraint is a kinematic (1) or driving (0) constraint
        myBodyI; % First body in constraint
        myBodyJ; % Second body in constraint
        myaBarI; % Vector defined on body I in body I reference frame
        myaBarJ; % Vector defined on body J in body J reference frame
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
        function obj = DP1constraint(bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot, constraintName)
            % Store name sent in by user
            if nargin < 8
                obj.myConstraintName = 'DP1 Constraint';
            else
                obj.myConstraintName = constraintName;
            end
            
            % Store attributes.
            obj.myBodyI = bodyI;
            obj.myBodyJ = bodyJ;
            obj.myaBarI = aBarI;
            obj.myaBarJ = aBarJ;
            obj.myFt = ft;
            obj.myFtDot = ftDot;
            obj.myFtDDot = ftDDot;
            
        end
        function obj = computeDP1constraint(obj, sys, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
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
                obj.computePhi(sys);
            end
            if(nuFlag == 1)
                obj.computeNu();
            end
            if(gammaFlag == 1)
                obj.computeGamma(sys);
            end
            if(phiPartialRFlag == 1)
                obj.computePhiPartialR(sys);
            end
            if(phiPartialPFlag == 1)
                obj.computePhiPartialP(sys);
            end
        end
        function obj = computePhi(obj,sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            aBarI = obj.myaBarI;
            aBarJ = obj.myaBarJ;
            ft = obj.myFt;
            t = obj.myTime;
            
            % Compute necessary parameters
            % Orientation matrix
            sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            sys.myBodies{bodyJ}.computeA();
            Aj = sys.myBodies{bodyJ}.myA;
            
            % Compute phi
            ftVal = ft(t);
            phi = aBarI'*Ai'*Aj*aBarJ - ftVal;
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
            aBarJ = obj.myaBarJ;
            t = obj.myTime;
            ftDDot = obj.myFtDDot;
            
            % Compute necessary parameters
            % Orientation matrix
            sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            sys.myBodies{bodyJ}.computeA();
            Aj = sys.myBodies{bodyJ}.myA;
            
            % Time derivatives of Euler parameters
            pDotI = sys.myBodies{bodyI}.myPDot;
            pDotJ = sys.myBodies{bodyJ}.myPDot;
            
            % B matrix
            sys.myBodies{bodyI}.computeB(aBarI);
            BmatrixI = sys.myBodies{bodyI}.myB;
            
            sys.myBodies{bodyJ}.computeB(aBarJ);
            BmatrixJ = sys.myBodies{bodyJ}.myB;
            
            % Time derivative of B matrix
            sys.myBodies{bodyI}.computeBDot(aBarI);
            BdotI = sys.myBodies{bodyI}.myBDot;
            
            sys.myBodies{bodyJ}.computeBDot(aBarJ);
            BdotJ = sys.myBodies{bodyJ}.myBDot;
            
            % a and aDot
            ai = Ai*aBarI;
            aDotI = BmatrixI*pDotI;
            
            aj = Aj*aBarJ;
            aDotJ = BmatrixJ*pDotJ;

            %Compute right hand side of acceleration equation
            gamma = -ai'*BdotJ*pDotJ - aj'*BdotI*pDotI - 2*aDotI'*aDotJ + ftDDot(t);
            obj.myGamma = gamma;
        end
        function obj = computePhiPartialR(obj, sys)
            % Extract both bodies
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            if (isGroundJ == 1)
                phiPartialRI = [0 0 0];
                phiPartialR = phiPartialRI;
            elseif (isGroundI == 1)
                phiPartialRJ = [0 0 0];
                phiPartialR = phiPartialRJ;
            else
                phiPartialRI = [0 0 0];
                phiPartialRJ = [0 0 0];
                phiPartialR = [phiPartialRI phiPartialRJ];
            end
            obj.myPhiPartialR = phiPartialR;
        end
        function obj = computePhiPartialP(obj,sys)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            aBarI = obj.myaBarI;
            aBarJ = obj.myaBarJ;
            
            % Check if either body is the ground
            isGroundI = sys.myBodies{bodyI}.myIsGround;
            isGroundJ = sys.myBodies{bodyJ}.myIsGround;
            
            % If bodyJ is the ground, this body contributes
            % nothing to the Jacobian
            if (isGroundJ == 1);
                % Compute orientation matrix, A, for bodyJ
                sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
                
                % Compute B(p, aBar) for bodyI
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                % Compute a for bodyI
                aj = Aj*aBarJ;
                
                % Compute phiPartialP
                phiPartialPI = aj'*BmatrixI;
                phiPartialP = phiPartialPI;
                
                % Special case when body I is ground
            elseif (isGroundI == 1);
                % Compute orientation matrix, A, for bodyI
                Ai = eye(3,3);
                
                % Compute B(p, aBar) for bodyI
                sys.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                
                % Compute a for bodyI
                ai = Ai*aBarI;
                
                % Compute phiPartialP
                phiPartialPJ = ai'*BmatrixJ;
                phiPartialP = phiPartialPJ;
            else
                % Compute orientation matrix, A, for each body
                sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                
                sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
                
                % Compute B(p, aBar) for bodyI and body J
                sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myB;
                
                sys.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = sys.myBodies{bodyJ}.myB;
                
                % Compute a for both bodies
                ai = Ai*aBarI;
                aj = Aj*aBarJ;
                
                % Compute phiPartialP
                phiPartialPI = aj'*BmatrixI;
                phiPartialPJ = ai'*BmatrixJ;
                phiPartialP = [phiPartialPI phiPartialPJ];
            end
            
            obj.myPhiPartialP = phiPartialP;
        end
    end
    
end
