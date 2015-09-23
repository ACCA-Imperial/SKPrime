classdef G03con < skpUnitTest.baseTestG0
% Test G0 in 3-connected domain.

% Copyright Everett Kropf, 2015
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
    D = skpDomain([-0.2517+0.3129i, 0.2307-0.4667i], ...
        [0.2377, 0.1557])
end

methods(Test)
    function simpleTest1(test)
        test.alpha = 0.5;
        checkViaSKProd(test)
    end
    
    function simpleTest2(test)
        test.alpha = -0.34 - 0.34i;
        checkViaSKProd(test)
    end
    
    function nearBdry1(test)
        [d, q] = domainData(test.D);
        test.alpha = d(1) + (q(1) + 1e-6);
        checkViaSKProd(test)
    end
    
    function alphaAtZero3con(test)
        test.alpha = 0;
        checkViaSKProd(test);
    end
end

end
