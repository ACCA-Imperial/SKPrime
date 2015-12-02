classdef schwarz < skpObject
%schwarz attempts to solve the modified Schwarz problem
%
% phi = schwarz(r, D, N)
%   Attempts to solve the Schwarz problem on multiply connected domain D
%   given the imaginary part of a function r(z) = imag(f) analytic on D.
%   The resulting object evaluates the solution of the Schwarz problem on
%   the boundary of the circles C_j. That is, for z_j on C_j, 
%
%     f(z_j) = phi(z_j) + 1i*r(z_j) + 1i*gamma_j
%
%   where the constant gamma_j is found as part of the solution.
%
% The current solution phi is given by the series
%
%                       -- N
%                       \
%   phi(z_j) = a(0,j) +  .     a(n,j)*eta(z_j)^n ...
%                       /                  + conj(a(n,j))*eta(z_j)^(-n)
%                       -- n=1
% where
%
%   eta(z_j) = (z_j - d_j)/q_j
%
% and a(0,0) = 0 since the BVP is determined up to a real constant. Note
% the constant gamma_j is the imaginary part of a(0,j).
%
% See also schwarzMatrix

% E. Kropf, 2015
% 
% This file is part of SKPrime.
% 
% SKPrime is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SKPrime is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SKPrime.  If not, see <http://www.gnu.org/licenses/>.

properties
    nTrapezoid = 100    % number of points for trapezoidal rule
end

properties(SetAccess=protected)
    theMatrix           % schwarzMatrix object
    phiCoef             % solution series coefficients
    rhsFun              % given function
    truncation          % series truncation level
end

properties(Dependent)
    rhsValue            % linear system RHS vector
end

properties(Access=protected)
    kNumber             % number of transform "modes" to use
    numBaseRows         % matrix construction value
end

methods
    function phi = schwarz(domain, N, nt)
        if ~nargin
            return
        end

        if isa(domain, 'schwarz')
            domain = domain.domain;
            if nargin < 2 && ~isempty(N)
                N = domain.truncation;
            end
            phi.theMatrix = domain.theMatrix;
        elseif isa(domain, 'circleRegion')
            domain = skpDomain(domain);
        end
        
        phi.domain = domain;
        phi.truncation = N;
        if nargin > 2
            phi.nTrapezoid = nt;
        end
        
        [~, ~, m] = domainDataB(domain);
        phi.kNumber = N + 1;
        phi.numBaseRows = (m + 1)*phi.kNumber - 1;
        
        if isempty(phi.theMatrix)
        	L = makeL(phi);
            Q = phi.numBaseRows;
            phi.theMatrix = ...
                schwarzMatrix([L; conj(L(:,Q+1:end)) conj(L(:,1:Q))]);
        end
    end
    
    function val = feval(phi, z)
        %provides function evaluation
        %
        % val = feval(phi, z)
        
        if isempty(phi.phiCoef)
            schwarz.noSolverWarn();
            val = nan(size(z));
            return
        end
        
        [d, q] = domainDataB(phi.domain);
        a = phi.phiCoef;
        N = (size(a, 1) - 1)/2;

        val = complex(nan(size(z)));
        for j = 1:size(a, 2)
            zj = z - d(j);
            L = abs(q(j) - abs(zj)) < eps(2);
            if any(L(:))
                val(L) = a(1,j) + ...
                    2*real(polyval([a(N+1:-1:2,j); 0], zj(L)/q(j)));
            end
        end
    end
    
    function phi = solve(phi, rfun)
        %call the solver for the modified Schwarz problem
        %
        % phi = solve(phi, rfun)
        %   Note this creates a new object.
        
        phi.rhsFun = rfun;
        phi.phiCoef = solveSystem(phi);
    end
    
    function r = get.rhsValue(phi)
        r = calcRHS(phi);
    end
end

methods(Access=protected)
    function b = calcRHS(phi)
        %computes the RHS of the linear system
        
        [d, q, m] = domainDataB(phi.domain);
        K = phi.kNumber;
        Q = phi.numBaseRows;
        rfun = phi.rhsFun;
        
        nt = phi.nTrapezoid;
        dt = 2*pi/nt;
        eit = exp(1i*(0:nt-1)*dt);
        zj = complex(zeros(1, (m + 1)*nt));
        rje = zj;
        for j = 0:m
            col = j*nt+(1:nt);
            zj(col) = d(j+1) + q(j+1)*eit;
            rje(col) = ((1 - (j > 0)*2)*q(j+1)*dt)*rfun(zj(col)).*eit;
        end

        b = complex(zeros(Q, (m+1)*nt));
        b(1,:) = rje;
        for k = 2:K-1
            b(k,:) = b(k-1,:).*zj;
        end
        for p = 1:m
            rowbase = p*K - 1;
            zjp = (zj - d(p+1))/q(p+1);
            b(rowbase+1,:) = rje./zjp;
            for k = 2:K
                b(rowbase+k,:) = b(rowbase+k-1,:)./zjp;
            end
        end
        b = sum(b, 2);
    end
    
    function L = makeL(phi)
        %constructs the matrix for the linear system
        
        [d, q, m] = domainDataB(phi.domain);
        N = phi.truncation;
        Q = phi.numBaseRows;
        
        L = complex(zeros(Q, 2*Q));
        for p = 0:m
            dp = d(p+1);
            qp = q(p+1);
            rows = (p > 0)*(N + (p - 1)*(N + 1)) + (1:N+(p>0));
            
            for j = 0:m
                dj = d(j+1);
                qj = q(j+1);
                cols = (j > 0)*(N + (j - 1)*(N + 1)) + (1:N+(j>0));
                
                if p == 0
                    cols = Q + cols;
                    if j == 0
                        B = diag(2i*pi*ones(N, 1));
                    else
                        if dj == 0
                            B = [zeros(N, 1), diag(-2i*pi*q(j+1).^(1:N))];
                        else
                            B = zeros(N, N + 1);
                            B(:,2) = -2i*pi*qj*dj.^(0:N-1);
                            for n = 3:N+1
                                B(n-1:N,n) = B(n-2:N-1,n-1)*qj.*(n-2:N-1)'/(n-2);
                            end
                        end
                    end
                else
                    if j == 0
                        if dp == 0
                            B = [zeros(1,N); diag(2i*pi*qp.^(2:N+1))];
                        else
                            B = zeros(N+1, N);
                            B(1,:) = 2i*pi*qp*dp.^(1:N);
                            B(2,1) = 2i*pi*qp^2;
                            for n = 2:N
                                brow = (2:n+1)';
                                k = -(brow - 1);
                                B(brow,n) = B(brow-1,n-1)*qp*(-n)./k;
                            end
                        end
                    elseif j == p
                        B = diag(-2i*pi*qp*ones(N + 1, 1));
                    else
                        cols = Q + cols;
                        djp = d(j+1) - d(p+1);
                        B = zeros(N+1);
                        k = (1:N+1)';
                        B(:,2) = -2i*pi*qp.^k*qj.*djp.^(-k);
                        for n = 3:N+1
                            B(:,n) = B(:,n-1)*qj/djp.*((-k-n+3)./(n-2));
                        end
                    end
                end
                
                L(rows,cols) = B;
            end
        end
    end
    
    function a = solveSystem(phi)
        %solves the linear system to find the coefficients
        
        [~, ~, m] = domainDataB(phi.domain);
        N = phi.truncation;
        M = 1 + 2*N;
        
        b = calcRHS(phi);
        x = phi.theMatrix\[b; conj(b)];

        a = zeros(M, m+1);
        a(2:end,1) = [x(1:N); conj(x(N:-1:1))];
        for j = 2:m+1
            rowbase = N+(j-2)*(N+1);
            a(:,j) = [x(rowbase+(1:N+1)); conj(x(rowbase+(N+1:-1:2)))];
        end
    end
end

methods(Access=private,Static)
    function noSolverWarn()
        warning('SKPrime:noSolverCalled', ...
                ['The solver has not been called for this "schwarz" ' ...
                'object. The function returned nan(size(z)).'])
    end
end

end
