classdef(Abstract) primeNewBase < skpUnitTest.skpTestBase
%primeNewBase is abstract base class for testing prime functions.

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

properties
    pointMap = containers.Map(...
        {'inner boundary', 'inner point',  'outer boundary', ...
         'outer point'}, ...
        {'innerBdryPoints', 'innerPoint', 'outerBdryPoints', ...
         'outerPoint'}, ...
        'UniformValues', true)
    
    innerBdryPoints
    innerPoint
    
    primeObj
    
    wprod
    prodLevel = 6
end

properties(Dependent)
    outerBdryPoints
    outerPoint
end

methods % getters
    function pts = get.outerBdryPoints(obj)
        pts = 1./conj(obj.innerBdryPoints);
    end
    
    function pt = get.outerPoint(obj)
        pt = 1/conj(obj.innerPoint);
    end
end

methods(TestClassSetup)
    function initTestPoints(test)
        test.innerPoint = test.domainData.testPointInside;
        test.innerBdryPoints = boundaryPts(test.domain, 5);
    end
    
    function createProduct(test)
        test.wprod = skprod(test.domain.dv, test.domain.qv, test.prodLevel);
    end
end

methods(TestMethodSetup)
    function createObjectForTest(test)
        test.primeObj = skprime(test.alpha, test.domain);
    end
end

methods(Test)
    function functionCheck(test)
        test.compareAllPoints(...
            @(z) test.wprod(z, test.alpha), test.primeObj, ...
            1e-3)
    end
end

end
