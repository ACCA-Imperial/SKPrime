classdef greensCjDa < bvpFun
%greensCjDa

% E. Kropf, 2016
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
    function da = greensCjDa(alpha, j, D, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 3
                sargs = {D, N};
            else
                sargs = {D};
            end
        end
        da = da@bvpFun(sargs{:});
        if ~nargin
            return
        end
  
        alpha = skpParameter(alpha, da.domain);
        da.parameter = alpha;

        % Known boundary part.
        da.singPart = @(z) 1i./(z - alpha)/(2*pi);
        h = @(z) -imag(da.singPart(z));
        
        % Solve.
        da.phiFun = solve(da.phiFun, h);
        
        % Normalisation.
        normval = da.phiFun.phiCoef(1,j+1);
        
        % Boundary function and domain interpolant.
        da.bdryFun = @(z) da.phiFun(z) - normval + 1i*h(z);
        da.contFun = SKP.bmcCauchy(da.bdryFun, da.domain, ...
            2*da.truncation);
    end
    
    function v = feval(da, z)
        v = complex(nan(size(z)));
        
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = da.hat(z(idx)) + da.singPart(z(idx));
        end
    end
    
    function v = hat(da, z)
        v = complex(nan(size(z)));
        z = z(:);
        
        % Points not "on" boundary.
        inD = true(size(z));
        inD(abs(1 - abs(z)) < eps(2)) = false;
        [d, q, m] = domainData(da.domain);
        for p = 1:m
            inD(abs(q(p) - abs(z - d(p))) < eps(2)) = false;
        end
        
        % Evaluate boundary points.
        if any(~inD(:))
            v(~inD) = da.bdryFun(z(~inD));
        end
        
        % Points not on boundary.
        if any(inD(:))
            v(inD) = da.contFun(z(inD));
        end
    end
end

end
