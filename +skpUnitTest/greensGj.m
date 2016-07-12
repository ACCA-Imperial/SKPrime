classdef(Abstract) greensGj < skpUnitTest.greensTestBase
%skpUnitTest.greensGj is the base class for testing Green's functions when
%j > 0.

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

properties(Abstract)
    jCircle
end

methods(TestMethodSetup)
    function createObjectForTest(test)
        test.gjObject = greensCj(test.alpha, test.jCircle, test.domain);
    end
    
    function specifyPotential(test)
        wp = test.wprod;
        j = test.jCircle;
        dj = test.domain.dv(j);
        qj = test.domain.qv(j);
        thj = @(z) test.domain.theta(j, z);
        
%         logprat = @(z,a) log(wp(z, a)./wp(z, thj(1/conj(a))))/2i/pi;
%         test.gjProd = @(z,a) logprat(z, a) + log(qj/abs(a - dj))/2i/pi;
        test.gjProd = @(z,a) log(wp(z, a)./wp(z, thj(1/conj(a))) ...
            *qj/abs(a - dj))/2i/pi;
        test.gjHatProd = @(z,a) ...
            test.gjProd(z, a) - log((z - a)./(z - thj(1/conj(a))))/2i/pi;
    end
end

methods(Test)
    % These override the usual checks until diff() and diffh() are added to
    % the greensCj class. Delete these overrides at that time.
    
    function hatVariableDerivative(test)
        gj = test.gjObject;
        try
            diffh(gj, 1)
            test.verifyFail(...
                'Calling ''diffh'' did not fail. Allow unit test.')
        catch
        end
    end
    
    function functionVariableDerivative(test)
        gj = test.gjObject;
        try
            diff(gj, 1)
            test.verifyFail(...
                'Calling ''diff'' did not fail. Allow unit test.')
        catch
        end
    end
end

end
