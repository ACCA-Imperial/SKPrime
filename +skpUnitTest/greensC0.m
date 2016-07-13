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
    
    gjRefFun
    gjRefHat
    
    gjProd
    gjHatProd
end

methods(TestClassSetup)
    function specifyReferenceProducts(test)
        wp = test.wprod;
        
        function g0 = g0prod(z, a)
            g0 = log(wp(z, a)./wp(z, 1/conj(a)))/2i/pi;
            if a ~= 0
                g0 = g0 - log(abs(a))/2i/pi;
            end
        end
        
        function g0h = g0hprod(z, a)
            g0h = g0prod(z, a);
            if a ~= 0
                g0h = g0h - log((z - a)./(z - 1/conj(a)))/2i/pi;
            else
                g0h = g0h - log(z)/2i/pi;
            end
        end

        test.gjProd = @g0prod;
        test.gjHatProd = @g0hprod;
    end
end

methods(TestMethodSetup)
    function createObjectForTest(test)
        test.gjObject = greensC0(test.alpha, test.domain);
    end
    
    function specifyReferenceFunctions(test)
        alpha = test.alpha;
        wt = skprime(alpha, test.domain);
        wb = invParam(wt);
        
        function g0 = g0fun(z)
            g0 = log(wt(z)./wb(z))/2i/pi;
            if alpha ~= 0
                g0 = g0 - log(abs(alpha))/2i/pi;
            end
        end
        
        function g0h = g0hat(z)
            g0h = g0fun(z);
            if alpha ~= 0
                g0h = g0h - log((z - alpha)./(z - 1/conj(alpha)))/2i/pi;
            else
                g0h = g0h - log(z)/2i/pi;
            end
        end
        
        test.gjRefFun = @g0fun;
        test.gjRefHat = @g0hat;
    end
end

end
