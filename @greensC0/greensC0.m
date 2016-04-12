classdef greensC0 < bvpFun
%greensC0 is the modified Green's function for MC domains wrt C0.
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
    normConstant = 0
    bdryFun
    contFun
end

methods
    function g0 = greensC0(alpha, D, N)
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
        aj = isclose(g0.domain, alpha);
        
        if ~isempty(alpha.ison) && alpha.ison == 0
            % Alpha on the unit circle means const zero.
            g0.logaFun = @(z) zeros(size(z));
            g0.singCorrFact = sf;
            g0.bdryFun = @(z) zeros(size(z));
            g0.contFun = @(z) zeros(size(z));
            return
        elseif alpha == 0
            for j = aj
                thj0 = d(j);
                thjinf = d(j) - q(j)^2/conj(d(j));
                sf = @(z) sf(z).*(z - thj0)./(z - thjinf);
            end            
            loga = @(z) log(z.*sf(z))/(2i*pi);
        elseif isinf(alpha)
            for j = aj
                thj0 = d(j);
                thjinf = d(j) - q(j)^2/conj(d(j));
                sf = @(z) sf(z).*(z - thjinf)./(z - thj0);
            end
            loga = @(z) log(sf(z)./z)/(2i*pi);
        else
            for j = isclose(g0.domain, alpha)
                thj = @(z) d(j) + q(j)^2*z./(1 - conj(d(j))*z);
                sf = @(z) sf(z).*(z - thj(alpha))./(z - thj(inv(alpha)));
            end
            
            loga = @(z) ...
                log((z - alpha)./(z - inv(alpha)).*sf(z))/(2i*pi);
        end
        g0.logaFun = loga;
        g0.singCorrFact = sf;
        
        if g0.domain.m == 0
            % Nothing more to do.
            return
        end
        
        % Known part on the boundary.
        ha = @(z) -imag(loga(z));
        
        % Solve for unknown part on the boundary.
        g0.phiFun = solve(g0.phiFun, ha);
        
        % Boundary function and Cauchy interpolant.
        g0.bdryFun = @(z) g0.phiFun(z) + 1i*ha(z);
        g0.contFun = SKP.bmcCauchy(g0.bdryFun, g0.domain, 2*g0.truncation);
        
        % Normalization factor.
        g0.normConstant = real(g0.g0hatEval(alpha) + log(sf(alpha))/(2i*pi));
    end
    
    function dg0 = diff(g0)
        %gives the derivative of the Green's function.
        %
        % dvj = diff(g0)
        %   Returns function handle to derivative of g0 function by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk.
        
        dg0h = diffh(g0);
        
        function dval = deval(z)
            alpha = g0.parameter;
            if ~isempty(alpha.ison) && alpha.ison == 0
                dval = complex(zeros(size(z)));
                return
            end
            
            dval = dg0h(z);
            if alpha == 0
                dval = dval + 1./z./(2i*pi);
            elseif isinf(alpha)
                dval = dval - 1./z./(2i*pi);
            else
                dval = dval + (1./(z - alpha) - 1./(z - conj(alpha)))/(2i*pi);
            end
        end
        
        dg0 = @deval;
    end
    
    function dgh = diffh(g0)
        %gives derivative of the analytic part wrt zeta variable.
        %
        % dvh = diffh(vj)
        %   Returns function handle to derivative of vj.hat by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk.
        
        dgh = dftDerivative(g0, @g0.hat);
    end
    
    function dgp = diffp(g0)
        %gives derivative with respect to parameter.
        %
        % dgp = diffp(g0)
        %   Derivative of g0 with respect to the parameter. The returned
        %   object has the parameter fixed and is a function of zeta.
        
        dgp = greensCjDa(g0.parameter, 0, g0);
    end
    
    function hatFun = G0hat(g0)
        hatFun = @g0.g0hatEval;
    end
    
    function v = g0hatEval(g0, z)
        %evaluate the "analytic" part.
        
        if g0.domain.m == 0
            v = complex(zeros(size(z)));
            return
        end
        
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
        v = complex(nan(size(z)));
        
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
        v = g0.logaFun(z) + g0.g0hatEval(z) - g0.normConstant;
    end
end

end
