classdef greensC0Dp < bvpFun
%greensC0Dp is the derivative with respect to the parameter of G0.
%
%  dpg0 = greensC0Dp(g0)
%    Solves the boundary value problems for the derivative with respect to
%    the parameter of the Greens function with respect to C0.
%
%  dpg0 = greensC0Dp(parameter, D)
%
%See also greensC0, bvpFun.

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
    
    partialWrtX
    partialWrtY
    normalizeConstant = 0
end

methods
    function dpg0 = greensC0Dp(g0, varargin)
        alpha = [];
        if ~nargin
            args = {};
        elseif (isa(g0, 'double') || isa(g0, 'skpParameter')) ...
                && nargin == 2 && (isa(varargin{1}, 'skpDomain') ...
                || isa(varargin{1}, 'bvpFun'))
            alpha = g0;
            args = varargin(1);
        elseif isa(g0, 'greensC0')
            alpha = g0.parameter;
            args = {g0};
        else
            error(PoTk.ErrorIdString.InvalidArgument, ...
                'Arguments not recognized.')
        end
        
        dpg0 = dpg0@bvpFun(args{:});
        if ~nargin
            return
        end
        
        alpha = skpParameter(alpha, dpg0.domain);
        dpg0.parameter = alpha;
        
        resx = @(z) (1./(alpha - z) + 1./(conj(alpha) ...
            - z*conj(alpha)^2))/2i/pi;
        dpg0.partialWrtX = ...
            genericPlusSingular(resx, @(z) -imag(resx(z)), dpg0);        
        
        resy = @(z) (1./(alpha - z) - 1./(conj(alpha) ...
            - z*conj(alpha)^2))/2/pi;
        dpg0.partialWrtY = ...
            genericPlusSingular(resy, @(z) -imag(resy(z)), dpg0);
        
        dpg0.normalizeConstant = dpg0.hat(alpha) + 1/(4i*pi*alpha);
    end
    
    function v = feval(dp, z)
        v = (dp.partialWrtX(z) - 1i*dp.partialWrtY(z))/2 ...
            - dp.normalizeConstant;
    end
    
    function v = hat(dp, z)
        %evalutate the "analytic" part of the function.
        
        v = (dp.partialWrtX.hat(z) - 1i*dp.partialWrtY.hat(z))/2 ...
            - dp.normalizeConstant;
    end
end

methods(Access=protected)
end

end

















