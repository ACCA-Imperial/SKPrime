classdef primeBasic < skpUnitTest.primeBase
%primeBasic unit tests on 3-connected domain.
%
% Tests: see "methods(Test)" section of source file.
%
% See also skpUnitTest.primeBase

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
    D = skpDomain([-0.2517+0.3129i, 0.2307-0.4667i], [0.2377, 0.1557])
end

methods(Test)
    function parameterInsideDisk(test)
        checkViaSKProd(test, -0.3 - 0.3i)
    end
    
    function parameterOutsideDisk(test)
        checkViaSKProd(test, 1/conj(-0.3 - 0.3i))
    end
    
    function parameterAtOrigin(test)
        checkViaSKProd(test, 0)
    end
    
    function parameterAtInfinity(test)
        checkViaSKProd(test, inf)
    end
    
    function parameterOnUnitCircle(test)
        checkViaSKProd(test, exp(2i*pi*rand(1)))
    end
    
    function parameterOnInnerBoundary(test)
        d = test.D.dv(1);
        q = test.D.qv(1);
        alpha = d + q*exp(5i*pi/4);
        checkViaSKProd(test, alpha)
    end
    
    function parameterOnOuterBoundary(test)
        d = test.D.dv(1);
        q = test.D.qv(1);
        alpha = 1/conj(d + q*exp(5i*pi/4));
        checkViaSKProd(test, alpha)
    end
    
    function parameterNearBoundary(test)
        d = test.D.dv(1);
        q = test.D.qv(1);
        checkViaSKProd(test, d + q + 1e-6)
    end
    
    function unitRootCheckInner(test)
        np = 200;
        zt = exp(2i*pi*(0:np-1)'/np);
        
        alpha1 = 0.6-0.192i;
        w1 = skprime(alpha1, test.D);
        alpha2 = 0.6-0.191i;
        w2 = skprime(alpha2, test.D);
        
        wt1 = w1(zt);
        wt2 = w2(zt);
        
        sameBranch = all(sign(real(wt1)) == sign(real(wt2))) ...
            & all(sign(imag(wt1)) == sign(imag(wt2)));
        verifyTrue(test, sameBranch, ...
            'Boundary values on different root branches.')
        verifyFalse(test, any(isnan(wt1) | isnan(wt2)), ...
            'Computation produced NaN values.')
    end
    
    function unitRootCheckOuter(test)
        np = 200;
        zt = exp(2i*pi*(0:np-1)'/np);
        zt = 1./conj(zt);
        
        alpha1 = 1/conj(0.6-0.192i);
        w1 = skprime(alpha1, test.D);
        alpha2 = 1/conj(0.6-0.191i);
        w2 = skprime(alpha2, test.D);
        
        wt1 = w1(zt);
        wt2 = w2(zt);
        
        sameBranch = all(sign(real(wt1)) == sign(real(wt2))) ...
            & all(sign(imag(wt1)) == sign(imag(wt2)));
        verifyTrue(test, sameBranch, ...
            'Boundary values on different root branches.')
        verifyFalse(test, any(isnan(wt1) | isnan(wt2)), ...
            'Computation produced NaN values.')
    end
end

end
