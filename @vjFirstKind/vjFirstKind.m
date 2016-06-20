classdef vjFirstKind < bvpFun
%vjFirstKind represents v_j function via Cauchy's theorem.
%
% vj = vjFirstKind(j, D)
% vj = vjFirstKind(j, D, N)
%   Solves the boundary value problem for the first-kind integral
%   associated with circle boundary 0 < j < m, where m is the number of
%   circles in the unit disk. The unit domain is specified by the skpDomain
%   object D. Truncation level of the Fourier series on each boundary is
%   given by N (see bvpFun.truncation for the default).
%
% See also: skpDomain, bvpFun

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

properties(SetAccess=protected)
    boundary
    
    logjFun
    bdryFun
    contFun
    
    thetaj
end

properties(Dependent)
    constants
    taujj
end

methods
    function vj = vjFirstKind(j, D, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 2
                sargs = {D, N};
            else
                sargs = {D};
            end
        end
        vj = vj@bvpFun(sargs{:});
        if ~nargin
            return
        end
        
        vj.boundary = j;
        [d, q] = domainData(vj.domain);
        vj.thetaj = @(z) d(j) + q(j)^2*z./(1 - conj(d(j))*z);
        
        if abs(d(j)) > q(j)
            % "double" circle center.
            dij = d(j)/(abs(d(j))^2 - q(j)^2);
            logj = @(z) log((z - d(j))./(z - dij))/(2i*pi);
        elseif abs(d(j)) > eps(2)
            logj = @(z) log(z - d(j))/(2i*pi);
        else
            logj = @(z) log(z)/(2i*pi);
        end
        vj.logjFun = logj;
        
        % Known part on the boundary.
        hj = @(z) -imag(logj(z));
        
        % Solve for unknown part on the boundary.
        vj.phiFun = solve(vj.phiFun, hj);
        
        % Boundary function and Cauchy interpolant.
        vj.bdryFun = @(z) vj.phiFun(z) + 1i*hj(z);
        vj.contFun = SKP.bmcCauchy(vj.bdryFun, vj.domain, vj.truncation);
    end
    
    function dvj = diff(vj)
        %gives first-kind integral derivative.
        %
        % dvj = diff(vj)
        %   Returns function handle to derivative of vj function by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk.
        
        dvh = diffh(vj);
        
        function dval = deval(z)
            j = vj.boundary;
            [dv, qv, ~, dvi] = domainData(vj.domain);
            
            dval = dvh(z) + 1./(z - dv(j))/(2i*pi);
            if abs(dv(j)) > qv(j)
                dval = dval - 1./(z - dvi(j))/(2i*pi);
            end
        end
        
        dvj = @deval;
    end
    
    function dvh = diffh(vj)
        %gives derivative of the analytic part.
        %
        % dvh = diffh(vj)
        %   Returns function handle to derivative of vj.hat by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk.
        
        dvh = dftDerivative(vj, @vj.hat);
    end
        
    function v = hat(vj, z)
        v = vjHatEval(vj, z);
    end
    
    function hatFun = vjHat(vj)
        hatFun = @vj.vjHatEval;
    end
    
    function v = vjHatEval(vj, z)
        v = zeros(size(z));
        z = z(:);
        
        % Points not "on" boundary.
        inD = true(size(z));
        inD(abs(1 - abs(z)) < eps(2)) = false;
        [d, q, m] = domainData(vj.domain);
        for p = 1:m
            inD(abs(q(p) - abs(z - d(p))) < eps(2)) = false;
        end
        
        % Evaluate boundary points.
        if any(~inD(:))
            v(~inD) = vj.bdryFun(z(~inD));
        end
        
        % Points not on boundary.
        if any(inD(:))
            v(inD) = vj.contFun(z(inD));
        end
    end
    
    function v = feval(vj, z)
        D = vj.domain;
        thj = vj.thetaj;
        tjj = vj.taujj;
        
        v = complex(nan(size(z)));
        notNan = ~isnan(z);
        
        % Points in D_zeta.
        mask = isin(D, z) & notNan;
        if any(mask(:))
            v(mask) = vj.logPlus(z(mask));
        end
        
        % Points in first reflection in to C_j.
        done = mask;
        mask(mask) = false;
        mask(~done) = isin(D, thj(1./conj(z(~done)))) & notNan(~done);
        if any(mask(:))
            v(mask) = conj(vj.logPlus(thj(1./conj(z(mask)))) - tjj);
        end
        
        % Points in D_zeta'.
        done = done | mask;
        mask(mask) = false;
        mask(~done) = isin(D, 1./conj(z(~done))) & notNan(~done);
        if any(mask(:))
            v(mask) = conj(vj.logPlus(1./conj(z(mask))));
        end
        
        % Points in first reflection into C_j'.
        done = done | mask;
        mask(mask) = false;
        mask(~done) = isin(D, thj(z(~done))) & notNan(~done);
        if any(mask(:))
            v(mask) = vj.logPlus(thj(z(mask))) - tjj;
        end
    end
    
    function v = logPlus(vj, z)
        v = vj.logjFun(z) + vj.vjHatEval(z);
    end
end

methods % Property access.
    function gjk = get.constants(vj)
        %Returns the constant values of imag(v_j) on circles C_k.
        
        gjk = imag(vj.phiFun.phiCoef(1,2:end));
    end
    
    function tau = get.taujj(vj)
        %retrieve the tau_jj constant associated with the function.
        
        tau = 2i*vj.constants(vj.boundary);
    end
end

end % vjFirstKind
