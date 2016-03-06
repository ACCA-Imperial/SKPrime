classdef skprimeDa < bvpFun
%skprimeDa

% E. Kropf, 2016
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

properties(SetAccess=protected)
    constant
    greensDeriv
    primeFun
end

methods
    function da = skprimeDa(w, N)
        if ~nargin
            sargs = {};
        else
            if nargin > 1
                sargs = {w.domain, N};
            else
                sargs = {w.domain};
            end
        end
        da = da@bvpFun(sargs{:});
        if ~nargin
            return
        end
        
        da.primeFun = w;
        alpha = w.parameter;
        dv = domainData(w.domain);
        
        [~, j] = max(abs([double(alpha); alpha - dv]));
        j = j - 1;
        assert(alpha ~= 0 || j > 0, ...
            'Something went horribly wrong!')
        
        if j == 0
            da.constant = 1./alpha;
        else
            da.constant = 1./(alpha - dv(j+1));
        end
        
        da.greensDeriv = greensCjDa(alpha, j, w);
        
        warning('SKP:warning', ...
            ['Prime parameter derivative is experimental. ', ...
            'Evaluation is acually the square of the derivative.'])
    end
    
    function v = feval(da, z)
        v = complex(nan(size(z)));
        
        inUnit = abs(z) <= 1 + eps(2);
        notNan = ~isnan(z);
        idx = inUnit & notNan;
        if any(idx(:))
            v(idx) = (4i*pi*da.greensDeriv(z(idx)) ...
                + da.constant).*da.primeFun.Xeval(z(idx));
        end
    end
end

end
