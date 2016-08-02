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

        if j < 1 || gj.domain.m < j
            error('The argument "j" must index an inner boundary circle.')
        end
        gj.boundary = j;
        
        alpha = skpParameter(alpha, gj.domain);
        gj.parameter = alpha;

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

        sf = SKP.SingularityCorrectionGj(alpha, j, gj.domain);
        gj.singCorrFact = sf;
        
        thja = gj.domain.theta(j, 1/conj(alpha));
        loga = @(z) log((z - alpha)./(z - thja).*sf(z))/2i/pi;
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
        gj.normConstant = real(gj.continuedFunction(alpha) ...
            + log(sf(alpha))/(2i*pi));
    end
    
    function dgj = diff(gj, n)
        %gives the variable derivative of the Green's function.
        %
        % dgj = diff(gj)
        % dgj = diff(gj, n)
        %   Returns function handle to derivative of g0 with respect to the
        %   variable by way of DFT on the boundary and Cauchy continuation
        %   for the interior. Derivative is restricted to the unit disk.
        %   The order of the derivative is given by n, which is computed
        %   by recursive application (default n=1).
        
        alpha = gj.parameter;
        if ~isempty(alpha.ison) && alpha.ison == gj.boundary
            dgj = @(z) complex(zeros(size(z)));
            return
        end
        
        if nargin < 2
            n = 1;
        end
        
        dgjh = diffh(gj, n);
        dmult = (-1)^(n-1)*factorial(n-1);
        dsf = diff(gj.singCorrFact);
        thja = gj.domain.theta(gj.boundary, 1/conj(alpha));

        dgj = @(z) dgjh(z) + ...
            (dmult*(1./(z - alpha).^n - 1./(z - thja).^n) + dsf(z))/2i/pi;
    end
    
    function dgh = diffh(gj, n)
        %gives derivative of the analytic part wrt zeta variable.
        %
        % dgh = diffh(gj)
        % dgh = diffh(gj, n)
        %   Returns function handle to derivative of g0.hat by way of
        %   DFT on the boundary and Cauchy continuation for the interior.
        %   Derivative is restricted to the unit disk. The order of the
        %   derivative is given by integer n > 0 (default = 1), computed by
        %   recursive application of the DFT.
        
        if nargin < 2
            n = 1;
        end
        
        dgh = dftDerivative(gj, @gj.hat, n);
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
        %provides function evaluation for the Green's function.
        
        v = gj.logaFun(z) + gj.hat(z);
    end
    
    function v = hat(gj, z)
        v = bvpEval(gj, z) - gj.normConstant;
    end
end

end
