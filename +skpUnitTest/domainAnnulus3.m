classdef domainAnnulus3 < skpUnitTest.domainForTest
%skpUnitTest.domainAnnulus3 is a 3-conneced annular domain.

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
    dv = [0; 0.52621+0.23899i]
    qv = [0.25; 0.14095]
    
    parameterInside = -0.50729+0.29388i
    testPointInside = 0.27638-0.55977i
end

methods
    function dom = domainAnnulus3()
        dom = dom@skpUnitTest.domainForTest();
        dom.parameterMap.remove('origin');
        dom.parameterMap.remove('infinity');
    end
end

methods(Static)
    function str = parameterLocations()
        str = skpUnitTest.domainForTest.parameterLocationsWithout(...
            'domainAnnulus3', 'origin', 'infinity');
    end
    
    function str = parameterLocationsWithout(varargin)
        str = skpUnitTest.domainForTest.parameterLocationsWithout(...
            'domainAnnulus3', varargin{:});
    end
end

end
