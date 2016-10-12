classdef basicConstraint
    %basicConstraint.m defines the basic kinematics constraints
    
    properties
        myConstraintName; % String. Name of the constraint.
        myConstraintType; % String. Type of basic constraint. Currently only DP1 or CD
        myConstraintAttributes; % Structure containing attributes of your selected constraint.
                                % These will be different depending on your
                                % selected constraint.
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
            if 
            
            
            
        end
    end
    
end

