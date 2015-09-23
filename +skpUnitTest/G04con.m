classdef G04con < skpUnitTest.baseTestG0
% Check G0 in 4-connected domain.

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
    D = skpDomain(...
        [0.052448+0.36539i, -0.27972-0.12762i, 0.48252-0.28147i], ...
        [0.15197, 0.17955, 0.20956])
end

methods(Test)
    function simpleDomainCheck(test)
        test.alpha = -0.43 + 0.29i;
        checkViaSKProd(test, 1e-4)
    end
    
    function nearBdry3(test)
        [d, q] = domainData(test.D);
        test.alpha = d(3) + (q(3) + 1e-6)*exp(5i*pi/4);
        checkViaSKProd(test)
    end
end

end
