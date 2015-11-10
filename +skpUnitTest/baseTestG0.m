classdef baseTestG0 < matlab.unittest.TestCase
% Base test routines for G0.

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

properties(Abstract)
    D
end

properties
    alpha
    zt
    wprod
    level = 7;
end

methods
    function T = baseTestG0(prodLvl)
        if ~nargin
            return
        end
        
        T.level = prodLvl;
    end
end

methods(TestClassSetup)
    function productFunctions(test)
        [dv, qv] = domainData(test.D);
        test.wprod = skprod(dv, qv, test.level);
    end
    
    function testPoints(test)
        [d, q, m] = domainDataB(test.D);
        
        np = 10;
        eit = exp(2i*pi*(0:np-1)'/np);
        z = complex(zeros(np, m+1));
        z(:,1) = 0.9*eit;
        z(:,2:end) = bsxfun(@plus, d(2:end).', ...
            bsxfun(@times, 1.15*q(2:end).', eit));
        
        test.zt = z(isin(test.D, z));
    end
end

methods
    function checkViaSKProd(test, reltol)
        if nargin < 2 || isempty(reltol)
            reltol = 1e-5;
        end
        
        g0 = greensC0(test.alpha, test.D);
        refval = refProd(test);
        relerr = abs(errval(test, refval, g0))./abs(refval);
        verifyLessThan(test, max(relerr), reltol)
    end
    
    function val = refProd(test)
        if test.alpha == 0 || isinf(test.alpha)
            val = log(test.wprod(test.zt, test.alpha)...
                ./test.wprod(test.zt, 1/conj(test.alpha)))/(2i*pi);
        else
            val = log(test.wprod(test.zt, test.alpha)...
                ./test.wprod(test.zt, 1/conj(test.alpha))...
                /abs(test.alpha))/(2i*pi);
        end
    end
    
    function err = errval(test, refval, g0)
        err = refval - g0(test.zt);
        rerr = real(err);
        % KLUDGE: gets around slight branch cut differences in the methods.
        Lerr = abs(rerr) > 0.9;
        rerr(Lerr) = rerr(Lerr) - sign(rerr(Lerr));
        err = complex(rerr - mean(rerr(:)), imag(err));
    end
end

end
