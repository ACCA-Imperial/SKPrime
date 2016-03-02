classdef greensC0da < bvpFun
%greensC0da

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
    parameter
    
    singPart
    bdryFun
    contFun
end

methods
    function dg = greensC0da(alpha, D, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 2
                sargs = {D, N};
            else
                sargs = {D};
            end
        end
        dg = dg@bvpFun(sargs{:});
        if ~nargin
            return
        end
  
        alpha = skpParameter(alpha, dg.domain);
        dg.parameter = alpha;
        
        % Known boundary bit.
        if isinf(alpha)
            dg.singPart = @(z) zeros(size(z));
        else
            dg.singPart = @(z) -1./(z - alpha)/(2*pi);
        end
        
        h = @(z) real(dg.singPart(z));
        
        % Solve.
        dg.phiFun = solve(dg.phiFun, h);
        
        % Boundary function and continuation.
        dg.bdryFun = @(z) dg.phiFun(z) + 1i*h(z);
        dg.contFun = SKP.bmcCauchy(dg.bdryFun, dg.domain, 2*dg.truncation);
    end
    
    function v = feval(dg, z)
        v = complex(nan(size(z)));
        
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = dg.hat(z(idx)) + dg.singPart(z(idx));
        end
    end
    
    function v = hat(dg, z)
        v = complex(nan(size(z)));
        z = z(:);
        
        % Points not "on" boundary.
        inD = true(size(z));
        inD(abs(1 - abs(z)) < eps(2)) = false;
        [d, q, m] = domainData(dg.domain);
        for p = 1:m
            inD(abs(q(p) - abs(z - d(p))) < eps(2)) = false;
        end
        
        % Evaluate boundary points.
        if any(~inD(:))
            v(~inD) = dg.bdryFun(z(~inD));
        end
        
        % Points not on boundary.
        if any(inD(:))
            v(inD) = dg.contFun(z(inD));
        end
    end
end

end
