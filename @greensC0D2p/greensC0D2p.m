classdef greensC0D2p < bvpFun
%greensC0D2p is the 2nd order derivative with respect to the 
%parameter of G0.
%
%  d2pg0 = greensC0D2p(g0)
%    Solves the boundary value problems for the 2nd order derivative with 
%    respect to the parameter of the Greens function with respect to C0.
%
%  d2pg0 = greensC0D2p(parameter, D)
%
%See also greensC0, bvpFun.

% Rhodri Nelson, 2016
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
    
    partialWrtO1
    % O1 refers to operator 1 = d^2/dx^2 - d^2/dy^2
    partialWrtO2
    % O2 refers to operator 2 = d^2/dxdy
    normalizeConstant = 0
end

methods
    function d2pg0 = greensC0D2p(g0, varargin)
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
        
        d2pg0 = d2pg0@bvpFun(args{:});
        if ~nargin
            return
        end
        
        alpha = skpParameter(alpha, d2pg0.domain);
        d2pg0.parameter = alpha;
        
        res1 = @(z) -1./(1i*pi).*(1./(z - alpha).^2 ...
               + 1./(conj(alpha).^2.*(conj(alpha).*z-1)) ...
               + z./(conj(alpha).*(conj(alpha).*z-1).^2));
        d2pg0.partialWrtO1 = ...
            genericPlusSingular(res1, @(z) -imag(res1(z)), d2pg0);        
        
        res2 = @(z) -1./(2*pi).*(1./(z - alpha).^2 ...
               - 1./(conj(alpha).^2.*(conj(alpha).*z-1)) ...
               - z./(conj(alpha).*(conj(alpha).*z-1).^2));
        d2pg0.partialWrtO2 = ...
            genericPlusSingular(res2, @(z) -imag(res2(z)), d2pg0);
        
        ncp = diffh(greensC0Dp(alpha, d2pg0.domain));
        d2pg0.normalizeConstant = d2pg0.hat(alpha) + ncp(alpha) ...
                                  - 1/(4i*pi*alpha^2);
    end
    
    function v = feval(dp, z)
        v = (dp.partialWrtO1(z) - 2i*dp.partialWrtO2(z))/4 ...
            - dp.normalizeConstant;
    end
    
    function v = hat(dp, z)
        %evalutate the "analytic" part of the function.
        
        v = (dp.partialWrtO1.hat(z) - 2i*dp.partialWrtO2.hat(z))/4 ...
            - dp.normalizeConstant;
    end
end

end
