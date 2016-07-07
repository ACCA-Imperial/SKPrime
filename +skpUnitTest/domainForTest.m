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

properties(Constant)
    parameterOrigin = 0
    parameterInfinity = inf
    parameterUnit = exp(2i*pi*rand(1))

    defaultParameterKeys = {...
        'inside', 'outside', 'origin', 'infinity', 'unit', ...
        'innerCirc1', 'outerCirc1', 'nearCirc1', ...
        'innerCircm', 'outerCircm'}

    defaultTestPointKeys = {...
        'inside', 'outside'}
end

properties(Dependent)
    defaultParameterValues
    defaultTestPointValues

    parameterOutside
    parameterInnerCirc1
    parameterOuterCirc1
    parameterNearCirc1
    parameterInnerCircm
    parameterOuterCircm
    
    testPointOutside
end

methods % Getters
    function pv = get.defaultParameterValues(td)
        keys = td.defaultParameterKeys;
        pv = td.keyToValue('parameter', keys);
    end
    
    function tv = get.defaultTestPointValues(td)
        keys = td.defaultTestPointKeys;
        tv = td.keyToValue('testPoint', keys);
    end
    
    function a = get.parameterOutside(td)
        a = 1/conj(td.parameterInside);
    end
    
    function a = get.parameterInnerCirc1(td)
        d = td.dv(1);
        q = td.qv(1);
        a = d + q*exp(5i*pi/4);
    end
    
    function a = get.parameterOuterCirc1(td)
        a = 1/conj(td.parameterInnerCirc1);
    end
    
    function a = get.parameterNearCirc1(td)
        d = td.dv(1);
        q = td.qv(1);
        a = d + q + 1e-6;
    end
    
    function a = get.parameterInnerCircm(td)
        m = numel(td.dv);
        d = td.dv(m);
        q = td.qv(m);
        a = d + q*exp(5i*pi/4);
    end
    
    function a = get.parameterOuterCircm(td)
        a = 1/conj(td.parameterInnerCircm);
    end
    
    function tp = get.testPointOutside(td)
        tp = 1/conj(td.testPointInside);
    end
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
        castr = skpUnitTest.(class).defaultParameterKeys;
        mask = false(size(castr));
        for i = 1:numel(varargin)
            mask = mask | strcmp(varargin{i}, castr);
        end
        castr = castr(~mask);
    end
    
    function castr = testPointLocations()
        castr = skpUnitTest.domainForTest.testPointMap.keys;
    end

    function vv = keyToValue(prefix, keys)
        vv = cell(size(keys));
        for i = 1:numel(vv)
            vstr = keys{i};
            vstr(1) = upper(vstr(1));
            vv{i} = [prefix, vstr];
        end
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
