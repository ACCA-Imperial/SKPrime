classdef primeBase < matlab.unittest.TestCase
%Base class for prime unit tests.

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
    wprod
    level = 8
   
    dv
    qv
    dvi
    qvi
    
    zin
    zbInner
    zbOuter
end

properties(Access=private)
    diagStr
end

methods(TestClassSetup)
    function prodFunction(test)
        [d, q, ~, di, qi] = domainData(test.D);

        test.dv = d;
        test.qv = q;
        test.dvi = di;
        test.qvi = qi;
        test.wprod = skprod(d, q, test.level);
    end
    
    function testPoints(test)
        d = test.dv;
        q = test.qv;
        np = 10;
        eit = exp(2i*pi*(0:np-1)/np);
        
        ztmp = [0.9*eit.', ...
            bsxfun(@plus, d, bsxfun(@times, 1.15*q, eit)).'];
        test.zin = ztmp(isin(test.D, ztmp));
        test.zbInner = [eit.', ...
            bsxfun(@plus, d, bsxfun(@times, q, eit)).'];
        test.zbOuter = bsxfun(@plus, test.dvi, ...
            bsxfun(@times, test.qvi, eit)).';
        j = find(isnan(test.dvi));
        if ~isempty(j)
            t = angle(d(j) + q(j)*exp(1i*angle(d(j))));
            t = linspace(t - pi/2, t + pi/2, np);
            test.zbOuter(:,j) = 1./conj(bsxfun(@plus, d(j), ...
                bsxfun(@times, q(j), exp(1i*t)))).';
        end
    end
end

methods
    function checkViaSKProd(test, alpha, reltol, notouter)
        if nargin < 3 || isempty(reltol)
            reltol = 1e-6;
        end
        
        diagString(test, sprintf(...
            'Alpha at %s (modulus = %.6f, arg = %.6f*pi)', ...
            num2str(alpha), abs(alpha), angle(alpha)/pi))
        wfun = skprime(alpha, test.D);
        wref = @(z) test.wprod(z, alpha);

        checkPoints(test, wfun, wref, test.zbInner, reltol, ...
            'on inner boundary')
        if nargin < 4 || ~strcmp(notouter, 'noouter')
            checkPoints(test, wfun, wref, test.zbOuter, reltol, ...
                'on outer boundary')
        end
        checkPoints(test, wfun, wref, test.zin, reltol, ...
            'in domain')
    end
end

methods(Access=protected)
    function checkPoints(test, wfun, wref, ztest, reltol, ptloc)
        wval = wfun(ztest(:));
        refval = wref(ztest(:));
        relerr = abs(refval - wval)./abs(refval);
        verifyLessThan(test, max(relerr), reltol, ...
            sprintf('Relative error tolerance exceeded %s.%s', ...
            ptloc, test.diagStr))
        verifyFalse(test, any(isnan(wval(:))), sprintf(...
            'Computation produced NaN values %s.%s', ...
            ptloc, test.diagStr))
    end
    
    function diagString(test, str)
        test.diagStr = sprintf('%s\n%s', test.diagStr, str);
    end
end

end
