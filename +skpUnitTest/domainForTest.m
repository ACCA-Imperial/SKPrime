classdef(Abstract) domainForTest
%domainForTest base test domain abstract class.

properties(Abstract)
    dv
    qv
    
    parameterInside
end

properties(Dependent)
    parameterOutside
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
                a = td.parameter.Outside;
                
            otherwise
                a = nan;
        end
    end
end

methods %getters
    function a = get.parameterOutside(td)
        a = 1/conj(td.parameterInside);
    end
end

end
