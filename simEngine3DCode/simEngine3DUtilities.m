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
        function [L U p q] = lucp(A,tol,pm_opt)
            %LUCP     LU factorization with complete pivoting.
            %
            % To compute the LU factorization under default settings:
            %
            %  [L U p q] = lucp(A)
            %
            % This produces a factorization such that L*U = A(p,q).  Vectors p and q
            % permute the rows and columns, respectively.
            %
            % The pivot tolerance can be controlled:
            %
            %  [L U p q] = lucp(A,tol)
            %
            % The algorithm will terminate if the absolute value of the pivot is less
            % than tol.
            %
            % Permutation matrices can be generated:
            %
            %  [L U P Q] = lucp(A,tol,'matrix')
            %  [L U P Q] = lucp(A,tol,'sparse')
            %
            % The first will generate full permutation matrices P and Q such that
            % L*U = P*A*Q.  The second generates sparse P and Q.
            %
            % If A is sparse, L and U will be sparse.
            %
            % This function works on non-square matrices.
            %
            % Input:
            %  A = matrix
            %  tol = pivot tolerance (defualt is 1e-10)
            %  pm_opt = permutation output options
            %         = 'vector' for permutation vectors, L*U = A(p,q), defualt
            %         = 'matrix' for full permutation matrices, L*U = P*A*Q
            %         = 'sparse' for sparse permutation matrices, L*U = P*A*Q
            %
            % Output:
            %  L = lower triangular factor
            %  U = upper triangular factor
            %  p = row permutation
            %  q = column permutation
            %
            % Reference:
            %  Algorithm 3.4.2, Matrix Computations, Golub and Van Loan. (3rd ed)
            %
            % Other Implementations:
            %  Gaussian Elimination using Complete Pivoting.  Alvaro Moraes.
            %  http://www.mathworks.com/matlabcentral/fileexchange/25758
            %
            %  Gauss elimination with complete pivoting.  Nickolas Cheilakos.
            %  http://www.mathworks.com/matlabcentral/fileexchange/13451
            %  (Does not work with rectangular A)
            %
            %  Rank Revealing Code.  Leslie Foster.
            %  http://www.math.sjsu.edu/~foster/rankrevealingcode.html
            %  (Uses mex libraries for computation)
            %
            
            
            %
            % 2010-03-28 (nwh) first version.
            % 2010-04-14 (nwh) added more references.
            % 2010-04-24 (nwh) perform final column swap so U is well conditioned
            %
            
            %
            % License: see license.txt.
            %
            
            % handle optional inputs
            if nargin < 2 || isempty(tol)
                tol = 1e-10;
            end
            
            if nargin < 3 || isempty(pm_opt)
                pm_opt = 'vector';
            end
            
            if strcmp(pm_opt,'vector')
                pm_flag = false;
                sp_flag = false;
            elseif strcmp(pm_opt,'matrix')
                pm_flag = true;
                sp_flag = false;
            elseif strcmp(pm_opt,'sparse')
                pm_flag = true;
                sp_flag = true;
            else
                error('lucp:invalidOption','''%s'' is an un recognized option.',pm_opt)
            end
            
            [n m] = size(A);
            
            % pivot vectors
            p = (1:n)';
            q = (1:m)';
            
            % temp storage
            rt = zeros(m,1); % row temp
            ct = zeros(n,1); % col temp
            t = 0; % scalar temp
            
            for k = 1:min(n-1,m)
                % determine pivot
                [cv ri] = max(abs(A(k:n,k:m)));
                [rv ci] = max(cv);
                rp = ri(ci)+k-1;
                cp = ci+k-1;
                
                % swap row
                t = p(k);
                p(k) = p(rp);
                p(rp) = t;
                rt = A(k,:);
                A(k,:) = A(rp,:);
                A(rp,:) = rt;
                
                % swap col
                t = q(k);
                q(k) = q(cp);
                q(cp) = t;
                ct = A(:,k);
                A(:,k) = A(:,cp);
                A(:,cp) = ct;
                
                if abs(A(k,k)) >= tol
                    rows = (k+1):n;
                    cols = (k+1):m;
                    A(rows,k) = A(rows,k)/A(k,k);
                    A(rows,cols) = A(rows,cols)-A(rows,k)*A(k,cols);
                else
                    % stop factorization if pivot is too small
                    break
                end
            end
            
            % final column swap if m > n
            if m > n
                % determine col pivot
                [cv ci] = max(abs(A(n,n:m)));
                cp = ci+n-1;
                
                % swap col
                t = q(n);
                q(n) = q(cp);
                q(cp) = t;
                ct = A(:,n);
                A(:,n) = A(:,cp);
                A(:,cp) = ct;
            end
            
            % produce L and U matrices
            % these are sparse if L and U are sparse
            l = min(n,m);
            L = tril(A(1:n,1:l),-1) + speye(n,l);
            U = triu(A(1:l,1:m));
            
            % produce sparse permutation matrices if desired
            if pm_flag
                p = sparse(1:n,p,1);
                q = sparse(q,1:m,1);
            end
            
            % produce full permutation matrices if desired
            if ~sp_flag
                p = full(p);
                q = full(q);
            end
        end
    end
end

