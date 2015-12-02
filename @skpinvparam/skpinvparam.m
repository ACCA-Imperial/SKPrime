classdef skpinvparam < skprime
%SKPINVPARAM is the inverted parameter prime function.
%
% w = skprime(...)
% winv = skpinvparam(w)
%
% Given the S-K prime function w(zeta, alpha), this class represents the
% inverted (w.r.t. the unit circle) parameter prime function
%
%   w(zeta, 1/conj(alpha)) = -zeta/conj(alpha)*conj(w(1/conj(zeta), alpha))
%
% for 0 < |alpha| < inf, and
%
%   w(zeta, 1/conj(alpha)) = conj(w(1/conj(zeta), alpha))
%
% otherwise.

% Everett Kropf, 2015
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
    function skpc = skpinvparam(skp)
        if ~nargin
            return
        end
        
        if ~isa(skp, 'skprime')
            error('SKPrime:invalidArgument', 'Expected a "skprime" object.')
        end
        
        skps = copyProperties(skp);
        pnames = fieldnames(skps);
        for k = 1:numel(pnames)
            skpc.(pnames{k}) = skps.(pnames{k});
        end
    end
    
    function Xc = Xeval(skpc, z)
        %gives the square of the prime function
        
        Xc = conj(Xeval@skprime(skpc, 1./conj(z)));
        alpha = skpc.parameter;
        if 0 ~= alpha && ~isinf(alpha)
            Xc = (z/conj(alpha)).^2.*Xc;
        else
            Xc = z.^2*Xc;
        end
    end
    
    function wc = feval(skpc, z)
        %provides function evaluation
        
        wc = conj(feval@skprime(skpc, 1./conj(z)));
        alpha = skpc.parameter;
        if 0 ~= alpha && ~isinf(alpha)
            wc = -z/conj(alpha).*wc;
        else
            wc = z.*wc;
        end
    end
    
    function alphac = parameter(skpc)
        %access inverted parameter
        
        alphac = inv(skpc.parameter);
    end
end

end
