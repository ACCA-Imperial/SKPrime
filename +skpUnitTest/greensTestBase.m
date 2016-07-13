classdef(Abstract) greensTestBase < skpUnitTest.skpTestBase
%skpUnitTest.greensTestBase is the base class for testing Green's
%functions.

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
        {'inner point', 'inner boundary'}, ...
        {'innerPoint', 'innerBdryPoints'}, ...
        'UniformValues', true)
    
    innerPoint
    innerBdryPoints

    wprod
    prodLevel = 6
end

properties(Abstract)
    gjObject
    
    gjRefRun
    gjRefHat
    
    gjProd
    gjHatProd
end

methods(TestClassSetup)
    function createProduct(test)
        test.wprod = skprod(test.domain.dv, test.domain.qv, test.prodLevel);
    end
    
    function initTestPoints(test)
        test.innerBdryPoints = boundaryPts(test.domain, 5);
        test.innerPoint = test.domainData.testPointInside;
    end
end

methods(Test)
    function hatCheck(test)
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.gjRefHat(z)), ...
            @(z) exp(2i*pi*test.gjObject.hat(z)), ...
            1e-6)
    end
    
    function functionCheck(test)        
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.gjRefFun(z)), ...
            @(z) exp(2i*pi*test.gjObject(z)), ...
            1e-6)
    end
    
    function hatVariableDerivative(test)
        gj = test.gjObject;
        
        dgjh = diffh(gj, 1);
        d2gjh = diffh(gj, 2);
        
        h = 1e-6;
        d2ref = @(z) (dgjh(z + h) - dgjh(z - h))/2/h;
        
        test.compareAllPoints(d2ref, d2gjh, 1e-4)
    end
    
    function functionVariableDerivative(test)
        gj = test.gjObject;
        
        dgj = diff(gj, 1);
        d2gj = diff(gj, 2);
        
        h = 1e-6;
        d2ref = @(z) (dgj(z + h) - dgj(z - h))/2/h;
        
        test.compareAllPoints(d2ref, d2gj, 1e-4)
    end
end

end
