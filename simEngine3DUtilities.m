classdef simEngine3DUtilities
    %simEngine3DUtilities.m 
    % This class provides a collection of functions that are commonly used
    % within simEngine3D.
    methods (Static)
        function [ matrix ] = skewSym( vector )
            x = vector(1);
            y = vector(2);
            z = vector(3);
            
            matrix = [0 -z y;
                z 0 -x;
                -y x 0];
        end
        function A = p2A(p) 
            % Orientation matrix A from euler parameters p
            %
            % Function input:
            %   p = [e0;e1;e2;e3] euler parameters
            %
            % output:
            %   A = [3x3] orientation matrix
            e0 = p(1);
            e1 = p(2);
            e2 = p(3);
            e3 = p(4);
            A = 2*[e0^2+e1^2-0.5, e1*e2-e0*e3, e1*e3+e0*e2;
                e1*e2+e0*e3, e0^2+e2^2-0.5, e2*e3-e0*e1;
                e1*e3-e0*e2, e2*e3+e0*e1, e0^2+e3^2-0.5];
        end
        function dij = computeDij(sys, bodyI, bodyJ, sBarIP, sBarJQ)
            % Computes a vector (dij) from point sBarIP to sBarJQ in the
            % global reference frame. 
            %
            % Function inputs:
            % sys : multibodySystem class
            %   Class representing a multibody system as defined using the
            %   mutlibodySystem class
            %
            % bodyI : int
            %   Integer representing the body that the tail of the vector
            %   dij lies on.
            %
            % bodyJ : int
            %   Integer representing the body that the head of the vector
            %   dij lies on.
            %
            % sBarIP : 3x1 double
            %   Vector defining the location of the tail of the vector dij
            %   in the body i reference frame
            % 
            % sBarJQ : 3x1 double
            %   Vector defining the location of the head of the vector dij
            %   in the body j reference frame  
            %
            % Function outputs:
            % dij : 3x1 double
            %   Vector from sBarIP to sBarJQ in the global reference frame.
            %   
            
            % Extract the current position and orientation matrix for each
            % body
            ri = sys.myBodies{bodyI}.myR;
            sys.myBodies{bodyI}.computeA;
            Ai = sys.myBodies{bodyI}.myA;
            
            rj = sys.myBodies{bodyJ}.myR;
            sys.myBodies{bodyJ}.computeA;
            Aj = sys.myBodies{bodyJ}.myA;
            
            % Compute dij
            dij = rj + Aj*sBarJQ - ri - Ai*sBarIP;
        end
        function dijDot = computeDijDot(sys, bodyI, bodyJ, sBarIP, sBarJQ)
            % Computes time derivative of a vector (dij) from point sBarIP to sBarJQ in the
            % global reference frame. 
            %
            % Function inputs:
            % sys : multibodySystem class
            %   Class representing a multibody system as defined using the
            %   mutlibodySystem class
            %
            % bodyI : int
            %   Integer representing the body that the tail of the vector
            %   dij lies on.
            %
            % bodyJ : int
            %   Integer representing the body that the head of the vector
            %   dij lies on.
            %
            % sBarIP : 3x1 double
            %   Vector defining the location of the tail of the vector dij
            %   in the body i reference frame
            % 
            % sBarJQ : 3x1 double
            %   Vector defining the location of the head of the vector dij
            %   in the body j reference frame  
            %
            % Function outputs:
            % dijDot : 3x1 double
            %   Time derivative of vector from sBarIP to sBarJQ in the 
            %   global reference frame.
            %   
            
            % Extract the current position time derivative for each body
            riDot = sys.myBodies{bodyI}.myRDot;
            rjDot = sys.myBodies{bodyJ}.myRDot;
            
            % Extract the current orientation time derivative for each body
            piDot = sys.myBodies{bodyI}.myPDot;
            pjDot = sys.myBodies{bodyJ}.myPDot;
            
            % Compute B-matrix for each body
            sys.myBodies{bodyI}.computeB(sBarIP);
            BmatrixI = sys.myBodies{bodyI}.myB;
            sys.myBodies{bodyJ}.computeB(sBarJQ);
            BmatrixJ = sys.myBodies{bodyJ}.myB;
            
            % Compute dijDot
            dijDot = rjDot + BmatrixJ*pjDot - riDot - BmatrixI*piDot;
        end
        
    end
    
end

