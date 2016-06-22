classdef greensC0Dpxy < genericPlusSingular
%greensC0Dpxy represents the partial derivative of G0 with respect to the
%parameter in the x or y direction.

% Everett Kropf, 2016
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
    xyIndicator
end

methods
    function dpg0 = greensC0Dpxy(alpha, xory, D)
        if ~nargin
            args = {};
        else
            if isa(D, 'skpDomain')
                dom = D;
            elseif isa(D, 'bvpFun')
                dom = D.domain;
            else
                error('SKPrime:InvalidParameter', ...
                    'Expedted a "skpDomain" object as third argument.')
            end
                
            alpha = skpParameter(alpha, dom);
            switch xory
                case 'x'
                    singPart = @(z) (1./(alpha - z) ...
                        + 1./(conj(alpha) - z*conj(alpha)^2))/2i/pi;
                    singDz = @(z) (1./(alpha - z).^2 ...
                        + 1./(1 - z*conj(alpha)).^2)/2i/pi;
                    
                case 'y'
                    singPart = @(z) (1./(alpha - z) ...
                        - 1./(conj(alpha) - z*conj(alpha)^2))/2/pi;
                    singDz = @(z) 1.0/(2*pi)*(1./(z-alpha).^2 ...
                        - 1./(conj(alpha).*(z-1./conj(alpha))).^2);
                    
                otherwise
                    error('SKPrime:InvalidParameter', ...
                        'Second argument must be ''x'' or ''y''.')
            end
            known = @(z) -imag(singPart(z));
            
            args = {singPart, known, D};
        end
        
        dpg0 = dpg0@genericPlusSingular(args{:});
        if ~nargin
            return
        end
        
        dpg0.singularPartDz = singDz;
        dpg0.parameter = alpha;
        dpg0.xyIndicator = xory;
    end
end

end
