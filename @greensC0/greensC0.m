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
end

methods
    function g0 = greensC0(alpha, D, N)
        if ~nargin
            sargs = {};
        else
            if abs(alpha) > 1 + eps(2)
                error('SKPrime:InvalidArgument', ...
                    'Parameter must have magnitude <= 1.')
            end
            
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
  
        g0.parameter = skpParameter(alpha, g0.domain);
        alpha = g0.parameter;
        [d, q] = domainData(g0.domain);
        
        % Singularity correction factor.
        sf = @(z) 1;
        aj = isclose(g0.domain, alpha);
        
        if ~isempty(alpha.ison) && alpha.ison == 0
            % Alpha on the unit circle means const zero.
            g0.logaFun = @(z) zeros(size(z));
            g0.singCorrFact = sf;
            g0.boundaryFunction = @(z) zeros(size(z));
            g0.continuedFunction = @(z) zeros(size(z));
            return
        elseif alpha == 0
            for j = aj
                thj0 = d(j);
                thjinf = d(j) - q(j)^2/conj(d(j));
                sf = @(z) sf(z).*(z - thj0)./(z - thjinf);
            end            
            loga = @(z) log(z.*sf(z))/(2i*pi);
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
        g0.boundaryFunction = @(z) g0.phiFun(z) + 1i*ha(z);
        g0.continuedFunction = SKP.bmcCauchy(...
            g0.boundaryFunction, g0.domain, 2*g0.truncation);
        
        % Normalization factor.
        g0.normConstant = real(g0.hat(alpha) + log(sf(alpha))/(2i*pi));
    end
    
    function dg0 = diff(g0, n)
        %gives the variable derivative of the Green's function.
        %
        % dg0 = diff(g0
        % dg0 = diff(g0, n)
        %   Returns function handle to derivative of g0 with respect to the
        %   variable by way of DFT on the boundary and Cauchy continuation
        %   for the interior. Derivative is restricted to the unit disk.
        %   The order of the derivative is given by n, which is computed
        %   by recursive application (default n=1).

        alpha = g0.parameter;
        if ~isempty(alpha.ison) && alpha.ison == 0
            dg0 = @(z) complex(zeros(size(z)));
            return
        end
            
        if nargin < 2
            n = 1;
        end
        
        dg0h = diffh(g0, n);
        dmult = (-1)^(n-1)*factorial(n-1);
        
        function dval = deval(z)
            dval = dg0h(z);
            if alpha == 0
                dval = dval + dmult./z.^n/2i/pi;
            else
                dval = dval + dmult*(1./(z - alpha).^n ...
                    - 1./(z - 1/conj(alpha)).^n)/2i/pi;
            end
        end
        
        dg0 = @deval;
    end
    
    function dgh = diffh(g0, n)
        %gives derivative of the analytic part wrt zeta variable.
        %
        % dgh = diffh(g0)
        % dgh = diffh(g0, n)
        %   Returns function handle to derivative of g0.hat by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk. The order of the
        %   derivative is given by integer n > 0 (default = 1), computed by
        %   recursive application of the DFT.
        
        if nargin < 2
            n = 1;
        end
        
        function dngh = rDftDiff(fun, n)
            if n > 1
                fun = rDftDiff(fun, n-1);
            end
            dngh = dftDerivative(g0, fun);
        end
        
        dgh = rDftDiff(@g0.hat, n);
    end
    
    function dgp = diffp(g0)
        %gives derivative with respect to parameter.
        %
        % dgp = diffp(g0)
        %   Derivative of g0 with respect to the parameter. The returned
        %   object has the parameter fixed and is a function of zeta.
        
        dgp = greensC0Dp(g0);
    end
    
    function v = feval(g0, z)
        %provides function evaluation for the Green's function.
        
        v = g0.logaFun(z) + g0.hat(z);
    end
    
    function v = hat(g0, z)
        %evaluate the "analytic" part of the function.
        %
        %  g0 = greensC0(...);
        %  v = g0.hat(z);
        
        if g0.domain.m == 0
            v = complex(zeros(size(z)));
            return
        end
        
        v = bvpEval(g0, z) - g0.normConstant;
    end
end

methods(Access=protected)
    function v = innerHat(g0, z)
        %hat function for points inside the unit disk.
        
        v = bvpEval(g0, z) - g0.normConstant;
    end
end

end
