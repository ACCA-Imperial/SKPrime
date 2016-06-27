classdef(Abstract) skpTestBase < matlab.unittest.TestCase
%skpUnitTest.skpTestBase is the abstract foundation for SKPrime unit tests.

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

properties(Abstract, MethodSetupParameter)
    parameterAt
end

properties(Abstract)
    domainData
    pointMap
end

properties
    domain
    alpha
end

methods(TestClassSetup)
    function domainSetup(test)
        testDomain = test.domainData;
        test.domain = skpDomain(testDomain);
    end
end

methods(TestMethodSetup)
    function determineParameter(test, parameterAt)
        test.alpha = test.domainData.parameter(parameterAt);
    end
end

methods
    function compareAllPoints(test, ref, fun, tol)
        keys = test.pointMap.keys();
        for i = 1:numel(keys)
            z = test.(test.pointMap(keys{i}));
            err = ref(z) - fun(z);
            test.verifyLessThan(max(abs(err(:))), tol, ...
                sprintf('Absolute error > %.1e on %s check.', tol, keys{i}))
        end
    end
end

end
