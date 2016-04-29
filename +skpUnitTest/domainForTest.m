classdef(Abstract) domainForTest
%domainForTest base test domain abstract class.

properties(Abstract)
    dv
    qv
    
    insideDisk
end

properties(Dependent)
    outsideDisk
end

methods
    function D = skpDomain(td)
        D = skpDomain(td.dv, td.qv);
    end
end

methods %getters
    function a = get.outsideDisk(td)
        a = 1/conj(td.insideDisk);
    end
end

end
