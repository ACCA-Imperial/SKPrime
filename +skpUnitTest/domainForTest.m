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

    defaultParameterKeys = {...
        'inside', 'outside', 'origin'}
    defaultParameterValues = {...
        'parameterInside', 'parameterOutside', 'parameterOrigin'}

    defaultTestPointKeys = {...
        'inside', 'outside'}
    defaultTestPointValues = {...
        'testPointInside', 'testPointOutside'}
end

properties(Access=protected)
    parameterMap
    testPointMap
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
end

methods(Static)
    function castr = parameterLocations()
        castr = skpUnitTest.domainForTest.defaultParameterKeys;
    end
    
    function castr = parameterLocationsWithout(class, varargin)
        if isempty(class)
            class = 'domainForTest';
        end
        castr = skpUnitTest.(class).parameterLocations;
        mask = false(size(castr));
        for i = 1:numel(varargin)
            mask = mask | strcmp(varargin{i}, castr);
        end
        castr = castr(~mask);
    end
    
    function castr = testPointLocations()
        castr = skpUnitTest.domainForTest.testPointMap.keys;
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
