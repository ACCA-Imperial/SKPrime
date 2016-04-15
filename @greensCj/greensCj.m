classdef greensCj < bvpFun
%GjCAUCHY is the modified Green's function for MC domains for Cj.
%
% The modified Green's function for multiply connected domains is given in
% terms of the Shottky-Klien prime function, w(z,a). That is
%
%                        /                                \
%                1       |    q_j           w(z,a)        |
%   G_j(z,a) = ----- log | --------- -------------------- |
%              2i*pi     | |a - d_j| w(z,th_j(1/conj(a))) |
%                        \                                /
%
% where w(z,a) is the S-K prime function, and th_j(z) is the reflection of
% a point through the unit circle, followed by a reflection through C_j.

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
    boundary
    
    logaFun
    singCorrFact
    normConstant
end

methods
    function gj = greensCj(alpha, j, D, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 3
                sargs = {D, N};
            else
                sargs = {D};
            end
        end
        gj = gj@bvpFun(sargs{:});
        if ~nargin
            return
        end
        
        gj.boundary = j;
        alpha = skpParameter(alpha, gj.domain);
        gj.parameter = alpha;
        [d, q] = domainData(gj.domain);
        
        if ~isempty(alpha.ison) && alpha.ison == j
            % Alpha on C_j is const zero.
            gj.logaFun = @(z) zeros(size(z));
            gj.boundaryFunction = @(z) zeros(size(z));
            gj.continuedFunction = @(z) zeros(size(z));
            return
        elseif isinf(alpha)
            error('SKPrime:invalidArgument', ...
                'The BVP for Gj with infinite alpha is not yet defined.')
        end
        
        sf = @(z) 1;
        ak = abs(alpha - d);
        ak = find(q + eps(2) < ak & ak < q + 0.1);
        for k = ak(ak ~= j)'
            thka = d(k) + q(k)^2/conj(alpha - d(k));
            sf = @(z) sf(z).*(z - d(k))./(z - thka);
        end
        if 0.9 < abs(alpha) && abs(alpha) < 1 - eps(2)
            sf = @(z) sf(z)./(z - 1/conj(alpha));
        end
        gj.singCorrFact = sf;
        
        thja = d(j) + q(j)^2/conj(alpha - d(j));
        loga = @(z) log((z - alpha)./(z - thja).*sf(z))/(2i*pi);
        gj.logaFun = loga;
        
        % Known part on the boundary.
        ha = @(z) -imag(loga(z));
        
        % Solve for unknown part on the boundary.
        gj.phiFun = solve(gj.phiFun, ha);
        
        % First normalization -- imag(Gj) is zero on C_j
        zeroCj = gj.phiFun.phiCoef(1,j+1);
        
        % Boundary function and Cauchy interpolant.
        gj.boundaryFunction = @(z) gj.phiFun(z) + 1i*ha(z) - zeroCj;
        gj.continuedFunction = SKP.bmcCauchy(...
            gj.boundaryFunction, gj.domain, 2*gj.truncation);
        
        % Second normalisztion -- adjusts the real part of the Schwarz
        % problem correctly.
        gj.normConstant = real(gj.continuedFunction(alpha) + log(sf(alpha))/(2i*pi));
    end
    
    function dgp = diffp(gj)
        %gives derivative with respect to parameter.
        %
        % dgp = diffp(g0)
        %   Derivative of g0 with respect to the parameter. The returned
        %   object has the parameter fixed and is a function of zeta.
        
        dgp = greensCjDa(gj.parameter, gj.boundary, gj);
    end
    
    function v = feval(gj, z)
        v = complex(nan(size(z)));
        
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = gj.logPlus(z(idx));
        end
        idx = ~inUnit & notNan;
        if any(idx(:))
            v(idx) = conj(gj.logPlus(1./conj(z(idx))));
        end
    end
    
    function v = hat(gj, z)
        v = bvpEval(gj, z);
    end
    
    function v = logPlus(gj, z)
        v = gj.logaFun(z) + gj.hat(z) - gj.normConstant;
    end
end

end
