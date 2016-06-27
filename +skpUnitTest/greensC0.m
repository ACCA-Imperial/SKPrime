classdef(Abstract) greensC0 < skpUnitTest.skpTestBase
%skpUnitTest.greensC0 is the test class for G0.

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

% properties(ClassSetupParameter)
%     domainInput = struct(...
%         'simple3', skpUnitTest.domainSimple3, ...
%         'annulus3', skpUnitTest.domainAnnulus3);
% end

properties
    pointMap = containers.Map(...
        {'inner point', 'inner boundary'}, ...
        {'innerPoint', 'innerBdryPoints'}, ...
        'UniformValues', true)
    
    innerPoint
    innerBdryPoints

    g0object
    
    wprod
    prodLevel = 6
    g0prod
    g0hatProd
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

methods(TestMethodSetup)
    function createObjectForTest(test)
        test.g0object = greensC0(test.alpha, test.domain);
    end
       
    function specifyPotential(test)
        wp = test.wprod;
        logprat = @(z,a) log(wp(z,a)./wp(z, 1/conj(a)))/2i/pi;
        if test.alpha ~= 0
            g0p = @(z,a) logprat(z, a) - log(abs(a))/2i/pi;
            g0hp = @(z,a) g0p(z, a) - log((z - a)./(z - 1/conj(a)))/2i/pi;
        else
            g0p = @(z,a) logprat(z, a);
            g0hp = @(z,a) g0p(z, a) - log(z)/2i/pi;
        end
        test.g0prod = g0p;
        test.g0hatProd = g0hp;
    end
end

methods(Test)
    function hatCheck(test)
        g0 = test.g0object;
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.g0hatProd(z, test.alpha)), ...
            @(z) exp(2i*pi*g0.hat(z)), 1e-4)
    end
    
    function functionCheck(test)
        g0 = test.g0object;
        
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.g0prod(z, test.alpha)), ...
            @(z) exp(2i*pi*g0(z)), 1e-4)
    end
    
    function hatVariableDerivative(test)
        g0 = test.g0object;
        
        dg0h = diffh(g0, 1);
        d2g0h = diffh(g0, 2);
        
        h = 1e-6;
        d2ref = @(z) (dg0h(z + h) - dg0h(z - h))/2/h;
        
        test.compareAllPoints(d2ref, d2g0h, 1e-4)
    end
    
    function functionVariableDerivative(test)
        g0 = test.g0object;
        
        dg0 = diff(g0, 1);
        d2g0 = diff(g0, 2);
        
        h = 1e-6;
        d2ref = @(z) (dg0(z + h) - dg0(z - h))/2/h;
        
        test.compareAllPoints(d2ref, d2g0, 1e-4)
    end
end

end
