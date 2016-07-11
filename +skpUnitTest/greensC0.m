classdef(Abstract) greensC0 < skpUnitTest.greensTestBase
%skpUnitTest.greensC0 specialises skpUnitTest.greensTestBase for the
%Green's function wrt C0.

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
    gjObject
    
    gjProd
    gjHatProd
end

methods(TestMethodSetup)
    function createObjectForTest(test)
        test.gjObject = greensC0(test.alpha, test.domain);
    end
       
    function specifyPotential(test)
        wp = test.wprod;
        logprat = @(z,a) log(wp(z,a)./wp(z, 1/conj(a)))/2i/pi;
        if test.alpha ~= 0
            g0p = @(z,a) logprat(z, a) - log(abs(a))/2i/pi;
            g0hp = @(z,a) g0p(z, a) - log((z - a)./(z - 1/conj(a)))/2i/pi;
        else
            g0p = @(z,a) logprat(z, a);
            g0hp = @(z,a) g0p(z, a) - log(z)/2i/pi;
        end
        test.gjProd = g0p;
        test.gjHatProd = g0hp;
    end
end

end
