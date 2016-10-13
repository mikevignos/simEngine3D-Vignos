classdef DP1constraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
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
            % allBodies : struct
            %   Data structure containing all of the bodies in the sys.
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
            if(nuFlag == 1)
                obj = obj.computeNu();
            end
            if(gammaFlag == 1)
                obj = obj.computeGamma(sys);
            end
            if(phiPartialRFlag == 1)
                obj = obj.computePhiPartialR();
            end
            if(phiPartialPFlag == 1)
                obj = obj.computePhiPartialP(sys);
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
            
            % Compute necessary parameters for body I
            % Orientation matrix
            sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            % Compute necessary parameters for body J. If body J is 0, this
            % is the ground. Hardcode in parameters for this case.
            if (bodyJ == 0)
                Aj = [1 0 0;
                    0 1 0;
                    0 0 1];
            else
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
            end
            
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
            
            % Compute necessary parameters for body I
            % Orientation matrix
            sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeA();
            Ai = sys.myBodies{bodyI}.myA;
            
            % Time derivatives of Euler parameters
            pDotI = sys.myBodies{bodyI}.myPDot;
            
            % B matrix
            sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeB(aBarI);
            BmatrixI = sys.myBodies{bodyI}.myBmatrix;
            
            % Time derivative of B matrix
            sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeBDot(aBarI);
            BdotI = sys.myBodies{bodyI}.myBDot;
            
            % a and aDot
            ai = Ai*aBarI;
            aDotI = BmatrixI*pDotI;
            
            % Compute necessary parameters for body J. Check for special
            % case when bodyJ = 0, which means bodyJ is the ground.
            if (bodyJ == 0)
                Aj = [1 0 0;
                    0 1 0;
                    0 0 1];
                pDotJ = zeros(4,1);
                BdotJ = zeros(3,4);
                aj = Aj*aBarJ;
                aDotJ = zeros(4,1); % This is all zeros because pDotJ is all zeros.
            else
                % Orientation matrix
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
                
                % Time derivatives of Euler parameters
                pDotJ = sys.myBodies{bodyJ}.myPDot;
                
                % Compute B(p, aBar)
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = sys.myBodies{bodyJ}.myBmatrix;
                
                % Compute B(pDot, aBar)
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeBDot(aBarJ);
                BdotJ = sys.myBodies{bodyJ}.myBDot;
                
                % Compute a and aDot for both bodies
                aj = Aj*aBarJ;
                aDotJ = BmatrixJ*pDotJ;
            end
            
            %Compute right hand side of acceleration equation
            gamma = -ai'*BdotJ*pDotJ - aj'*BdotI*pDotI - 2*aDotI'*aDotJ + ftDDot(t);
            obj.myGamma = gamma;
        end
        function obj = computePhiPartialR(obj)
            % If bodyJ is the ground (i.e. bodyJ = 0), this body contributes
            % nothing to the Jacobian
            bodyJ = obj.myBodyJ;
            if (bodyJ == 0)
                phiPartialR = [0 0 0];
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
            
            % If bodyJ is the ground (i.e. bodyJ = 0), this body contributes
            % nothing to the Jacobian
            if (bodyJ == 0);
                % Compute orientation matrix, A, for bodyJ
                Aj = eye(3,3);
                                
                % Compute B(p, aBar) for bodyI
                sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myBmatrix;
                
                % Compute a for bodyI
                aj = Aj*aBarJ;
                
                % Compute phiPartialP
                phiPartialPI = aj'*BmatrixI;
                phiPartialP = phiPartialPI;
                
            else
                % Compute orientation matrix, A, for each body
                sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeA();
                Ai = sys.myBodies{bodyI}.myA;
                
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeA();
                Aj = sys.myBodies{bodyJ}.myA;
                
                % Compute B(p, aBar) for bodyI and body J
                sys.myBodies{bodyI} = sys.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = sys.myBodies{bodyI}.myBmatrix;
                
                sys.myBodies{bodyJ} = sys.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = sys.myBodies{bodyJ}.myBmatrix;
                
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
