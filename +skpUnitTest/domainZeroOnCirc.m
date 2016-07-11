classdef domainZeroOnCirc < skpUnitTest.domainForTest
%skpUnitTest.domainZeroOnCirc provides a domain with the origin on a
%circle.

properties
    dv = [-0.2-0.2i; 0.56+0.36i]
    qv = [abs(-0.2-0.2i); 0.18]
    
    parameterInside = -0.11545 + 0.4688i
    testPointInside = 0.39534 - 0.28688i
end

methods
    function dom = domainZeroOnCirc()
        dom = dom@skpUnitTest.domainForTest();
        dom.parameterMap.remove('origin');
        dom.parameterMap.remove('infinity');
    end
end

methods(Static)
    function str = parameterLocations()
        str = skpUnitTest.domainOffset3.parameterLocationsWithout(...
            'domainZeroOnCirc', 'origin', 'infinity');
    end
    
    function str = parameterLocationsWithout(varargin)
        str = skpUnitTest.domainForTest.parameterLocationsWithout(...
            'domainZeroOnCirc', varargin{:});
    end
end

end
