classdef bodies
    %BODIES.m Class that defines the attributes of each body in a multibody
    %system
    %   Detailed explanation goes here
    
    properties
        myBodyName; % Name of the body
        myBodyType; % Currently only bars (links) are supported
        myMass; % Mass of the body
        myLength; % Length of the body. This will only be set if the body is a bar.
        myPoints; % Structure containing the important points on the body.
        myVectors; % Structure containing important vectors on the body
    end
    
    methods
        function obj = bodies(bodyName, bodyType, mass, length)
            % Set parameters that are sent in
            
            % Compute polar moment of inertia for body
            
            
        end
    end
    
end

