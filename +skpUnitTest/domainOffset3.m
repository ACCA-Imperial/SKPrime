classdef domainOffset3 < skpUnitTest.domainForTest
%skpUnitTest.domainOffset3 is an offset annular domain.

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
    dv = [-0.1-0.1i; 0.56+0.36i]
    qv = [0.3; 0.18]
    
    parameterInside = -0.3 + 0.3i
    testPointInside = -0.17143 + 0.54577i;
end

methods
    function dom = domainOffset3()
        dom = dom@skpUnitTest.domainForTest();
        dom.parameterMap.remove('origin');
        dom.parameterMap.remove('infinity');
    end
end

methods(Static)
    function str = parameterLocations()
        str = skpUnitTest.domainOffset3.parameterLocationsWithout(...
            'domainOffset3', 'origin', 'infinity');
    end
    
    function str = parameterLocationsWithout(varargin)
        str = skpUnitTest.domainForTest.parameterLocationsWithout(...
            'domainOffset3', varargin{:});
    end
end

end
