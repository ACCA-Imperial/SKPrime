classdef primeSimple3_np < matlab.unittest.TestCase
%skpUnitTest.primeSimple3_np has unparameterized test cases for simple 3
%domain.

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
    domainData = skpUnitTest.domainSimple3
end

methods
    function D = skpDomain(test)
        D = skpDomain(...
            test.domainData.dv, test.domainData.qv);
    end
end

methods(Test)
    function unitRootCheckInner(test)
        D = skpDomain(test);
        np = 200;
        zt = exp(2i*pi*(0:np-1)'/np);
        
        alpha1 = 0.6-0.192i;
        w1 = skprime(alpha1, D);
        alpha2 = 0.6-0.191i;
        w2 = skprime(alpha2, D);
        
        wt1 = w1(zt);
        wt2 = w2(zt);
        
        sameBranch = all(sign(real(wt1)) == sign(real(wt2))) ...
            & all(sign(imag(wt1)) == sign(imag(wt2)));
        verifyTrue(test, sameBranch, ...
            'Boundary values on different root branches.')
        verifyFalse(test, any(isnan(wt1) | isnan(wt2)), ...
            'Computation produced NaN values.')
    end
    
    function unitRootCheckOuter(test)
        D = skpDomain(test);
        np = 200;
        zt = exp(2i*pi*(0:np-1)'/np);
        zt = 1./conj(zt);
        
        alpha1 = 1/conj(0.6-0.192i);
        w1 = skprime(alpha1, D);
        alpha2 = 1/conj(0.6-0.191i);
        w2 = skprime(alpha2, D);
        
        wt1 = w1(zt);
        wt2 = w2(zt);
        
        sameBranch = all(sign(real(wt1)) == sign(real(wt2))) ...
            & all(sign(imag(wt1)) == sign(imag(wt2)));
        verifyTrue(test, sameBranch, ...
            'Boundary values on different root branches.')
        verifyFalse(test, any(isnan(wt1) | isnan(wt2)), ...
            'Computation produced NaN values.')
    end
end

end