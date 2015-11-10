classdef vjCauchy < bvpFun
%vjCauchy represents v_j function via Cauchy's theorem.

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
end

methods
    function vj = vjCauchy(j, D, N)
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
        v = nan(size(z));
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = vj.logPlus(z(idx));
        end
        idx = ~inUnit & notNan;
        if any(idx(:))
            v(idx) = conj(vj.logPlus(1./conj(z(idx))));
        end
    end
    
    function v = logPlus(vj, z)
        v = vj.logjFun(z) + vj.vjHatEval(z);
    end
end

end % vjCauchy
