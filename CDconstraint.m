<<<<<<< HEAD
classdef CDconstraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
        myBodyI;
        myBodyJ;
        myCoordVec;
        mysBarIP;
        mysBarJQ;
        myFt;
        myFtDot;
        myFtDDot;
        myTime; % Value of current time step
        
        % Properties of constraint that can be computed
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates
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
        function obj = computeCDconstraint(obj, system, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
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
                obj = obj.computePhi(system);
            end
            if (nuFlag == 1)
                obj = obj.computeNu();
            end
            if (gammaFlag == 1)
                obj = obj.computeGamma(system);
            end
            if (phiPartialRFlag == 1)
                obj = obj.computePhiPartialR();
            end
            if (phiPartialPFlag == 1)
                obj = obj.computePhiPartialP(system);
            end
        end
        function obj = computePhi(obj,system)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            ft = obj.myFt;
            t = obj.myTime;
            
            % Compute necessary parameters for body I
            % Orientation matrix
            system.myBodies{bodyI} = system.myBodies{bodyI}.computeA();
            Ai = system.myBodies{bodyI}.myA;
            ri = system.myBodies{bodyI}.myR;
            
            % Compute necessary parameters for body J. If body J is 0, this
            % is the ground. Hardcode in parameters for this case.
            if (bodyJ == 0)
                Aj = eye(3,3);
                rj = zeros(3,1);
            else
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeA();
                Aj = system.myBodies{bodyJ}.myA;
                rj = system.myBodies{bodyJ}.myR;
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
        function obj = computeGamma(obj, system)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            t = obj.myTime;
            ftDDot = obj.myFtDDot;
            
            % Compute necessary parameters for body I         
            % Time derivatives of Euler parameters
            pDotI = system.myBodies{bodyI}.myPDot;
            
            % Time derivative of B matrix
            system.myBodies{bodyI} = system.myBodies{bodyI}.computeBDot(sBarIP);
            BdotI = system.myBodies{bodyI}.myBDot;
                        
            % Compute necessary parameters for body J. Check for special
            % case when bodyJ = 0, which means bodyJ is the ground.
            if (bodyJ == 0)
                pDotJ = zeros(4,1);
                BdotJ = zeros(3,4);
            else              
                % Time derivatives of Euler parameters
                pDotJ = system.myBodies{bodyJ}.myPDot;
                
                % Compute B(pDot, aBar)
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeBDot(sBarJQ);
                BdotJ = system.myBodies{bodyJ}.myBDot;
            end
            
            %Compute right hand side of acceleration equation
            ftDDotVal = ftDDot(t);
            gamma = coordVec'*BdotI*pDotI - coordVec*BdotJ*pDotJ + ftDDotVal;
            obj.myGamma = gamma;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%% STOPPED HERE  %%%%%%%%%%%%%%%%%%%%%
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
        function obj = computePhiPartialP(obj,system)
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
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = system.myBodies{bodyI}.myBmatrix;
                
                % Compute a for bodyI
                aj = Aj*aBarJ;
                
                % Compute phiPartialP
                phiPartialPI = aj'*BmatrixI;
                phiPartialP = phiPartialPI;
                
            else
                % Compute orientation matrix, A, for each body
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeA();
                Ai = system.myBodies{bodyI}.myA;
                
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeA();
                Aj = system.myBodies{bodyJ}.myA;
                
                % Compute B(p, aBar) for bodyI and body J
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = system.myBodies{bodyI}.myBmatrix;
                
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = system.myBodies{bodyJ}.myBmatrix;
                
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

=======
classdef CDconstraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
        myBodyI;
        myBodyJ;
        myCoordVec;
        mysBarIP;
        mysBarJQ;
        myFt;
        myFtDot;
        myFtDDot;
        myTime; % Value of current time step
        
        % Properties of constraint that can be computed
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates
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
        function obj = computeCDconstraint(obj, system, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
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
                obj = obj.computePhi(system);
            end
            if (nuFlag == 1)
                obj = obj.computeNu();
            end
            if (gammaFlag == 1)
                obj = obj.computeGamma(system);
            end
            if (phiPartialRFlag == 1)
                obj = obj.computePhiPartialR();
            end
            if (phiPartialPFlag == 1)
                obj = obj.computePhiPartialP(system);
            end
        end
        function obj = computePhi(obj,system)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            ft = obj.myFt;
            t = obj.myTime;
            
            % Compute necessary parameters for body I
            % Orientation matrix
            system.myBodies{bodyI} = system.myBodies{bodyI}.computeA();
            Ai = system.myBodies{bodyI}.myA;
            ri = system.myBodies{bodyI}.myR;
            
            % Compute necessary parameters for body J. If body J is 0, this
            % is the ground. Hardcode in parameters for this case.
            if (bodyJ == 0)
                Aj = eye(3,3);
                rj = zeros(3,1);
            else
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeA();
                Aj = system.myBodies{bodyJ}.myA;
                rj = system.myBodies{bodyJ}.myR;
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
        function obj = computeGamma(obj, system)
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            sBarIP = obj.mysBarIP;
            sBarJQ = obj.mysBarJQ;
            coordVec = obj.myCoordVec;
            t = obj.myTime;
            ftDDot = obj.myFtDDot;
            
            % Compute necessary parameters for body I         
            % Time derivatives of Euler parameters
            pDotI = system.myBodies{bodyI}.myPDot;
            
            % Time derivative of B matrix
            system.myBodies{bodyI} = system.myBodies{bodyI}.computeBDot(sBarIP);
            BdotI = system.myBodies{bodyI}.myBDot;
                        
            % Compute necessary parameters for body J. Check for special
            % case when bodyJ = 0, which means bodyJ is the ground.
            if (bodyJ == 0)
                pDotJ = zeros(4,1);
                BdotJ = zeros(3,4);
            else              
                % Time derivatives of Euler parameters
                pDotJ = system.myBodies{bodyJ}.myPDot;
                
                % Compute B(pDot, aBar)
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeBDot(sBarJQ);
                BdotJ = system.myBodies{bodyJ}.myBDot;
            end
            
            %Compute right hand side of acceleration equation
            ftDDotVal = ftDDot(t);
            gamma = coordVec'*BdotI*pDotI - coordVec*BdotJ*pDotJ + ftDDotVal;
            obj.myGamma = gamma;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%% STOPPED HERE  %%%%%%%%%%%%%%%%%%%%%
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
        function obj = computePhiPartialP(obj,system)
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
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = system.myBodies{bodyI}.myBmatrix;
                
                % Compute a for bodyI
                aj = Aj*aBarJ;
                
                % Compute phiPartialP
                phiPartialPI = aj'*BmatrixI;
                phiPartialP = phiPartialPI;
                
            else
                % Compute orientation matrix, A, for each body
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeA();
                Ai = system.myBodies{bodyI}.myA;
                
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeA();
                Aj = system.myBodies{bodyJ}.myA;
                
                % Compute B(p, aBar) for bodyI and body J
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeB(aBarI);
                BmatrixI = system.myBodies{bodyI}.myBmatrix;
                
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeB(aBarJ);
                BmatrixJ = system.myBodies{bodyJ}.myBmatrix;
                
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

>>>>>>> 65497c4dd5107ecd5f81f147ff75c57183ddcbb4
