classdef primeZeroOnBdry < skpUnitTest.primeBase
%primeZeroOnBdry unit tests on 3-connected domain.
%
% Put the origin on a boundary.
%
% The test for the parameter at the origin is conspicuously absent, since
% the inverse of the parameter would be at infinity on a boundary, which
% is not a covered case. So we make no claim for any functionality in this
% case.

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
    D = skpDomain([-0.2-0.2i, 0.56+0.36i], [abs(-0.2-0.2i), 0.18])
end

methods(Test)
    function parameterInsideDisk(test)
        checkViaSKProd(test, 0.6-0.4i, 1e-4)
    end
    
    function parameterOutsideDisk(test)
        checkViaSKProd(test, 1/conj(0.6-0.4i), 1e-4)
    end
    
    function parameterNearOrigin(test)
        % Near origin and boundary.
        alpha = 0.09*exp(1i*angle(-test.D.dv(1)));
        checkViaSKProd(test, alpha, 3e-6)
    end
    
    function parameterNearOriginBoundary(test)
        % Near origin boundary, away from origin.
        d = test.D.dv(1);
        q = test.D.qv(1);
        alpha = d + (q + 0.09)*exp(1i*pi/2);
        checkViaSKProd(test, alpha, 1e-5)
    end
    
    function parameterOnUnitCircle(test)
        checkViaSKProd(test, exp(2i*pi*rand(1)), 4e-5)
    end
    
    function parameterOnOriginBoundary(test)
        d = test.D.dv(1);
        q = test.D.qv(1);
        checkViaSKProd(test, d + q*exp(5i*pi/4), 1e-5, 'noouter')
    end
    
    function parameterOnOtherBoundary(test)
        d = test.D.dv(2);
        q = test.D.qv(2);
        checkViaSKProd(test, d + q*exp(5i*pi/4), 2e-5, 'noouter')
    end
end

end
