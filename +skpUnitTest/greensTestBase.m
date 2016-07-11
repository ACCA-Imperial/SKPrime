classdef(Abstract) greensTestBase < skpUnitTest.skpTestBase
%skpUnitTest.greensTestBase is the base class for testing Green's
%functions.

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
        gj = test.gjObject;
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.gjHatProd(z, test.alpha)), ...
            @(z) exp(2i*pi*gj.hat(z)), 1e-4)
    end
    
    function functionCheck(test)
        gj = test.gjObject;
        
        test.compareAllPoints(...
            @(z) exp(2i*pi*test.gjProd(z, test.alpha)), ...
            @(z) exp(2i*pi*gj(z)), 1e-4)
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
