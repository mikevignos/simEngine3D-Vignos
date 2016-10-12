classdef basicConstraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType = 'DP1'; % String. Type of basic constraint. Currently only DP1 or CD
        myConstraintAttributes; % Structure containing attributes of your selected constraint.
        % These will be different depending on your
        % selected constraint.
        % List of attributes needed for each constraint
        % DP1 constraint:
        % bodyI: ADD IN DESCRIPTION OF ATTRIBUTES
        % bodyJ
        % aBarI
        % aBarJ
        % ft : Function prescribing the value of the constraint
        % ftDot : Function prescribing the derivative of the constraint
        % ftDDot : Function prescribing the second derivative of the constraint
        %
        % CD constraint:
        % bodyI
        % bodyJ
        % coordVec
        % sBarIP
        % sBarJQ
        % ft
        myTime; % Value of current time step
        myPhi; % Value of the expression of the constraint at current time step
        myNu; % Right hand side of the velocity equation at current time step
        myGamma; % Right hand side of the acceleration equation at the current time step
        myPhiPartialR; % Partial derivative of phi w.r.t. the location generalized coordinates
        myPhiPartialP; % Partial derivative of phi w.r.t. the orientation generalized coordinates
    end
    
    methods
        function obj = basicConstraint(constraintName, constraintType, attributes)
            % Store name and constraint type sent in by user
            obj.myConstraintName = constraintName;
            obj.myConstraintType = constraintType;
            
            % Store attributes.
            obj.myConstraintAttributes = attributes;
            
            % Throw error to user if all necessary attributes are not
            % defined for constraint.
            switch constraintType
                case 'DP1'
                    requiredAttributes = [{'bodyI'} {'bodyJ'} {'aBarI'} {'aBarJ'} {'ft'} {'ftDot'} {'ftDDot'}];
                    for iA = 1:length(requiredAttributes)
                        if ~isfield(attributes,requiredAttributes{iA})
                            disp(['ERROR: ' requiredAttributes{iA} ' must be defined for ' constraintType ' constraint'])
                        end
                    end
                    
                case 'CD'
                    requiredAttributes = [{'bodyI'} {'bodyJ'} {'coordVec'} {'sBarIP'} {'sBarJQ'} {'ft'}];
                    for iA = 1:length(requiredAttributes)
                        if ~isfield(attributes,requiredAttributes{iA})
                            disp(['ERROR: ' requiredAttributes{iA} ' must be defined for ' constraintType ' constraint'])
                        end
                    end
            end
        end
        
        function obj = computeDP1constraint(obj, allBodies, t, phiFlag, nuFlag, gammaFlag, phiPartialRFlag, phiPartialPFlag)
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
            bodyI = obj.myAttributes.bodyI;
            bodyJ = obj.myAttributes.bodyJ;
            aBarI = obj.myAttributes.aBarI;
            aBarJ = obj.myAttributes.aBarJ;
            ft = obj.myAttributes.ft;  
            ftDot = obj.myAttributes.ftDot;
            ftDDot = obj.myAttributes.ftDDot;
            
            % Compute orientation matrix, A, for each body
            allBodies{bodyI} = allBodies{bodyI}.computeA();
            Ai = allBodies{bodyI}.myA;
            
            allBodies{bodyJ} = allBodies{bodyJ}.computeA();
            Aj = allBodies{bodyJ}.myA;
           
            % Compute 
            
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
                
                % Compute B(pDot
                BdotJ = 
                BdotI = 
                pDot = 
                gamma = 
            end
            if(phiPartialRFlag == 1)
            end
            if(phiPartialPFlag == 1)
            end
            
            
            
        end
    end
    
end

