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

methods
    function dpg0 = greensC0Dp(g0, varargin)
        if ~nargin
            args = {};
        else
            args = {things};
        end
        
        dpg0 = dpg0@bvpFun(args{:});
        if ~nargin
            return
        end
        
        dpg0 = solveBVPs(dpg0);
    end
    
    function v = feval(dp, z)
        
    end
end

methods(Access=protected)
    function dp = solveBVPs(dp)
        resx = @(z) (1./(alpha - z) + 1./(conj(alpha) - z*conj(alpha)^2))/2i/pi;
        phix = solve(g0.phiFun, @(z) -imag(resx(z)));
        pxg0h = @(z) phix(z) - 1i*imag(resx(z));
        pxg0 = @(z) pxg0h(z) + resx(z);
        
        
        resy = @(z) (1./(alpha - z) - 1./(conj(alpha) - z*conj(alpha)^2))/2/pi;
        phiy = solve(g0.phiFun, @(z) -imag(resy(z)));
        pyg0h = @(z) phiy(z) - 1i*imag(resy(z));
        pyg0 = @(z) pyg0h(z) + resy(z);
    end
end

end

















