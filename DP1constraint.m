classdef DP1constraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
        myBodyI;
        myBodyJ;
        myaBarI;
        myaBarJ;
        myFt;
        myFtDot;
        myFtDDot;
        myTime; % Value of current time step
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates
    end
    
    methods
        function obj = DP1constraint(constraintName, bodyI, bodyJ, aBarI, aBarJ, ft, ftDot, ftDDot)
            % Store name sent in by user
            obj.myConstraintName = constraintName;
            
            % Store attributes.
            obj.myBodyI = bodyI;
            obj.myBodyJ = bodyJ;
            obj.myaBarI = aBarI;
            obj.myaBarJ = aBarJ;
            obj.myFt = ft;
            obj.myFtDot = ftDot;
            obj.myFtDDot = ftDDot;
            
        end
        
        function obj = computeDP1constraint(obj, system, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
            % Computes necessary quantities for DP1 constraint
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
            
            % Extract needed attributes
            bodyI = obj.myBodyI;
            bodyJ = obj.myBodyJ;
            aBarI = obj.myaBarI;
            aBarJ = obj.myaBarJ;
            ft = obj.myFt;
            ftDot = obj.myFtDot;
            ftDDot = obj.myFtDDot;
            
            % Compute orientation matrix, A, for each body
            system.myBodies{bodyI} = system.myBodies{bodyI}.computeA();
            Ai = system.myBodies{bodyI}.myA;
            
            system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeA();
            Aj = system.myBodies{bodyJ}.myA;
                        
            % Compute desired parameters
            if (phiFlag == 1)
                ftVal = ft(t);
                phi = aBarI'*Ai'*Aj*aBarJ - ftVal;
                obj.myPhi = phi;
            end
            if(nuFlag == 1)
                nu = ftDot(t);
                obj.myNu = nu;
            end
            if(gammaFlag == 1)
                ai = Ai*aBarI;
                aj = Aj*aBarJ;
                pDotI = system.myBodies{bodyI}.myPDot;
                pDotJ = system.myBodies{bodyJ}.myPDot;
                
                % Compute B(pDot, aBar) for bodyI and bodyJ
                system.myBodies{bodyI} = system.myBodies{bodyI}.computeBDot(aBarI);
                BdotI = system.myBodies{bodyI}.myBDot;
                
                system.myBodies{bodyJ} = system.myBodies{bodyJ}.computeBDot(aBarJ);
                BdotJ = system.myBodies{bodyJ}.myBDot;
                
                
                gamma =
            end
            if(phiPartialRFlag == 1)
            end
            if(phiPartialPFlag == 1)
            end
            
            
            
        end
    end
    
end

