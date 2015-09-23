classdef G0Cauchy < bvpFun
%G0CAUCHY is the modified Green's function for MC domains for C0.
%
% The modified Green's function for multiply connected domains is given in
% terms of the Shottky-Klien prime function, w(z,a). That is
%
%                        /                   \
%                1       |      w(z,a)       |
%   G_0(z,a) = ----- log | ----------------- |
%              2i*pi     | |a|w(z,1/conj(a)) |
%                        \                   /.

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
        
    logaFun
    singCorrFact
    bdryFun
    contFun
end

methods
    function g0 = G0Cauchy(alpha, D, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 2
                sargs = {D, N};
            else
                sargs = {D};
            end
        end
        g0 = g0@bvpFun(sargs{:});
        if ~nargin
            return
        end
  
        alpha = skpParameter(alpha, g0.domain);
        g0.parameter = alpha;
        [d, q] = domainData(g0.domain);
        
        % Singularity correction factor.
        sf = @(z) 1;
        if ~isempty(alpha.ison) && alpha.ison == 0
            % Alpha on the unit circle is const zero.
            g0.logaFun = @(z) zeros(size(z));
            g0.singCorrFact = sf;
            g0.bdryFun = @(z) zeros(size(z));
            g0.contFun = @(z) zeros(size(z));
            return
        elseif alpha == 0
            aj = find(abs(d) < q + 0.1);
            for j = aj'
                acj = d(j) + q(j)^2/conj(-d(j));
                sf = @(z) sf(z).*(z - d(j))./(z - acj);
            end
            
            loga = @(z) log(z.*sf(z))/(2i*pi);
        elseif isinf(alpha)
            loga = @(z) -log(z)/(2i*pi);
        else
            if abs(alpha) <= 1
                aj = abs(alpha - d);
                aj = find(eps(2) < aj & aj < q + 0.1);
                for j = aj'
                    acj = d(j) + q(j)^2/conj(alpha - d(j));
                    sf = @(z) sf(z).*(z - d(j))./(z - acj);
                end
            else            
                aj = find(abs(1/conj(alpha) - d) < q + 0.1);
                for j = aj'
                    acj = d(j) + q(j)^2*alpha/(1 - conj(d(j))*alpha);
                    sf = @(z) sf(z).*(z - acj)./(z - d(j));
                end
            end
            
            loga = @(z) ...
                log((z - alpha)./(z - 1/conj(alpha)).*sf(z))/(2i*pi);
        end
        g0.logaFun = loga;
        g0.singCorrFact = sf;
        
        % Known part on the boundary.
        ha = @(z) -imag(loga(z));
        
        % Solve for unknown part on the boundary.
        g0.phiFun = solve(g0.phiFun, ha);
        
        % Boundary function and Cauchy interpolant.
        g0.bdryFun = @(z) g0.phiFun(z) + 1i*ha(z);
        g0.contFun = bmcCauchy(g0.bdryFun, g0.domain, 2*g0.truncation);
    end
    
    function hatFun = G0hat(g0)
        hatFun = @g0.g0hatEval;
    end
    
    function v = g0hatEval(g0, z)
        v = complex(nan(size(z)));
        z = z(:);
        
        % Points not "on" boundary.
        inD = true(size(z));
        inD(abs(1 - abs(z)) < eps(2)) = false;
        [d, q, m] = domainData(g0.domain);
        for p = 1:m
            inD(abs(q(p) - abs(z - d(p))) < eps(2)) = false;
        end
        
        % Evaluate boundary points.
        if any(~inD(:))
            v(~inD) = g0.bdryFun(z(~inD));
        end
        
        % Points not on boundary.
        if any(inD(:))
            v(inD) = g0.contFun(z(inD));
        end
    end
    
    function v = feval(g0, z)
        v = nan(size(z));
        
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = g0.logPlus(z(idx));
        end
        idx = ~inUnit & notNan;
        if any(idx(:))
            v(idx) = conj(g0.logPlus(1./conj(z(idx))));
        end
    end
    
    function v = hat(g0, z)
        v = g0hatEval(g0, z);
    end
    
    function v = logPlus(g0, z)
        v = g0.logaFun(z) + g0.g0hatEval(z);
    end
end

end
