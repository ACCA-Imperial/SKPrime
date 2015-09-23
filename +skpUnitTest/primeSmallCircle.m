classdef primeSmallCircle < skpUnitTest.primeBase
%primeSmallCircle tests for small interior circles.

% Copyright Everett Kropf, 2015
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

properties
    D = skpDomain([-0.5, 0.5], [1e-4, 1e-4])
end

methods(Test)
    function pointsOnBoundary(test)
        [d, q] = domainDataB(test.D);
        np = 200;
        zt = bsxfun(@plus, d, ...
            bsxfun(@times, q, exp(2i*pi*(0:np-1)/np))).';
        
        alpha = 0.5i;
        w = skprime(alpha, test.D);
        
        verifyFalse(test, any(isnan(w(zt(:)))), ...
            'Calculation produced NaN values.')
    end
end

end
