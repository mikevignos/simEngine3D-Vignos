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
        function p = A2p(A)
            % Compute Euler parameters from current orientation matrix
            %
            % input:
            %   A = [3x3] orientation matrix
            % output:
            %   p = [e0;e1;e2;e3] euler parameters
            
            e0 = sqrt((trace(A)+1)/4); % the sign of e0 is arbitrary
            
            if e0 ~= 0
                e1 = (A(3,2)-A(2,3))/(4*e0);
                e2 = (A(1,3)-A(3,1))/(4*e0);
                e3 = (A(2,1)-A(1,2))/(4*e0);
                
            elseif e0 == 0
                disp('Not implemented yet, see slide 25 (9/21/16) or page 341')
                return
            end
            
            p =[e0;e1;e2;e3];
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
        function B = computeBmatrix(p, aBar)
            e0 = p(1);
            e = p(2:4);
            eTilde = simEngine3DUtilities.skewSym(e);
            aBarTilde = simEngine3DUtilities.skewSym(aBar);
            
            % Create B matrix
            B1 = (e0*eye(3,3) + eTilde)*aBar;
            B2 = e*aBar' - (e0*eye(3,3) + eTilde)*aBarTilde;
            B = 2*[B1 B2];
        end
        function Kmatrix = computeKmatrix(aBar, b)
            bTilde = simEngine3DUtilities.skewSym(b);
            aBarTilde = simEngine3DUtilities.skewSym(aBar);
            
            K1 = aBar'*b;
            K2 = aBar'*bTilde;
            K3 = aBarTilde*b;
            K4 = aBar*b' + b*aBar' - aBar'*b*[1 0 0; 0 1 0; 0 0 1];
            K = [K1 K2;
                K3 K4];
            
            Kmatrix = 2*K;
        end
        function [pDot] = omegaBar2pDot(sys, bodyNumber, omegaBar)
            % Convert angular velocity in the body reference frame
            % (omegaBar) to pDot.
            %
            % Function inputs:
            % sys : structure
            %   Multibody system that contains the body for which you are
            %   computing pDot.
            %
            % bodyNumber : int
            %   Number of the body in the multibody system for which you
            %   are converting omegaBar to pDot
            %
            % omegaBar : 3 x 1 vector
            %   Angular velocity in the body reference frame.
            %
            % Function outputs:
            % pDot : 4 x 1 vector
            %   First time derivative of Euler parameters for this body.
            
            body = bodyNumber;
            
            % Compute the Gmatrix for this body.
            sys.myBodies{body}.computeGmatrix();
            G = sys.myBodies{body}.myG;
            
            % Compute pDot
            pDot = 0.5*G'*omegaBar(:);
        end
%         function [x,U]=gausselim(A,b)
%             % function to perform gauss eliminination
%             %FORWARD ELIMINATION
%             n=length(b);
%             m=zeros(n,1);
%             x=zeros(n,1);
%             for k =1:n-1;
%                 %compute the kth column of M
%                 m(k+1:n) = A(k+1:n,k)/A(k,k);
%                 %compute
%                 %An=Mn*An-1;
%                 %bn=Mn*bn-1;
%                 for i=k+1:n
%                     A(i, k+1:n) = A(i,k+1:n)-m(i)*A(k,k+1:n);
%                 end;
%                 b(k+1:n)=b(k+1:n)-b(k)*m(k+1:n);
%             end
%             U= triu(A);
%             %BACKWARD ELIMINATION
%             x(n)=b(n)/A(n,n);
%             for k =n-1:-1:1;
%                 b(1:k)=b(1:k)-x(k+1)* U(1:k,k+1);
%                 x(k)=b(k)/U(k,k);
%             end
%         end
    end
end

