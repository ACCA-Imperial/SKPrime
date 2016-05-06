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

properties(Dependent)
    parameterOutside
    testPointOutside
end

properties(Constant)
    parameterOrigin = 0
end

properties(Access=protected)
    parameterMap
    defaultParameterKeys = {...
        'inside', 'outside', 'origin'}
    defaultParameterValues = {...
        'parameterInside', 'parameterOutside', 'parameterOrigin'}
    
    testPointMap
    defaultTestPointKeys = {...
        'inside', 'outside'}
    defaultTestPointValues = {...
        'testPointInside', 'testPointOutside'}
end

methods
    function dom = domainForTest()
        dom = useDefaultParameterMap(dom);
        dom = useDefaultTestPointMap(dom);
    end
    
    function D = skpDomain(td)
        D = skpDomain(td.dv, td.qv);
    end
    
    function a = get.parameterOutside(td)
        a = 1/conj(td.parameterInside);
    end
    
    function tp = get.testPointOutside(td)
        tp = 1/conj(td.testPointInside);
    end
    
    function out = subsref(td, S)
        if numel(S) == 2 && strcmp(S(1).type, '.') ...
                && numel(S(2).subs) == 1 && strcmp(S(2).type, '()')
            isat = S(2).subs{1};
            switch S(1).subs
                case 'parameter'
                    out = td.(td.parameterMap(isat));
                    return
                    
                case 'testPoint'
                    out = td.(td.testPointMap(isat));
                    return
            end
        end
        
        out = builtin('subsref', td, S);
    end
    
    function castr = parameterLocations(dom)
        castr = dom.parameterMap.keys;
    end
    
    function castr = testPointLocations(dom)
        castr = dom.testPointMap.keys;
    end
end

methods(Access=protected)
    function dom = useDefaultParameterMap(dom)
        dom.parameterMap = containers.Map(...
            dom.defaultParameterKeys, dom.defaultParameterValues);
    end
    
    function dom = useDefaultParameterMapMinusOrigin(dom)
        dom = useDefaultParameterMap(dom);
        remove(dom.parameterMap, 'origin');
    end
    
    function dom = useDefaultTestPointMap(dom)
        dom.testPointMap = containers.Map(...
            dom.defaultTestPointKeys, dom.defaultTestPointValues);
    end
end

end
