classdef(Abstract) domainForTest
%domainForTest base test domain abstract class.

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

properties(Abstract)
    dv
    qv
    
    parameterInside
    testPointInside
end

methods
    function D = skpDomain(td)
        D = skpDomain(td.dv, td.qv);
    end
    
    function a = parameter(td, isat)
        switch isat
            case 'inside'
                a = td.parameterInside;
                
            case 'outside'
                a = 1/conj(td.parameterInside);
                
            otherwise
                a = nan;
        end
    end
    
    function tp = testPoint(td, isat)
        switch isat
            case 'inside'
                tp = td.testPointInside;
                
            case 'outside'
                tp = 1/conj(td.testPointInside);
                
            otherwise
                tp = nan;
        end
    end
end

end
