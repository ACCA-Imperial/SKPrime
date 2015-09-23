classdef primeOffAnn < skpUnitTest.primeBase
%prime function tests on 3-connected domain.

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
    D = skpDomain([-0.1-0.1i, 0.56+0.36i], [0.3, 0.18])
end

methods(Test)
    function alphaSmallInFD(test)
        checkViaSKProd(test, -0.3 + 0.3i, 1e-4)
    end
    
    function alphaLargeInFD(test)
        checkViaSKProd(test, 1/conj(-0.3 + 0.3i), 1e-4)
    end
    
    function unitParameter(test)
        checkViaSKProd(test, exp(2i*pi*rand(1)), 9e-6)
    end
    
    function alphaOnAnnBdry(test)
        d = test.D.dv(1);
        q = test.D.qv(1);
        checkViaSKProd(test, d + q*exp(5i*pi/4), 1e-4)
    end
    
    function alphaOnOtherBdry(test)
        d = test.D.dv(2);
        q = test.D.qv(2);
        checkViaSKProd(test, d + q*exp(5i*pi/4), 1e-4)
    end
end

end
