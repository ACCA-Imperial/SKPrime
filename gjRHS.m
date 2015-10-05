classdef gjRHS < skpObject
%gjRHS is RHS with only Gj.

% E. Kropf, 2015
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
    parameter
    g0fun
    gjFuns
end

methods
    function rhs = gjRHS(param, D)
        if ~nargin
            return
        end
        
        if ~isa(param, 'skpParameter')
            error('SKPrime:invalidArgument', ...
                'Parameter must be a "skpParameter" object.')
        end
        rhs.parameter = param;
        rhs.domain = D;
        
        rhs.g0fun = G0Cauchy(param, D);
        rhs.gjFuns = cell(D.m, 1);
        for j = 1:D.m
            rhs.gjFuns{j} = GjCauchy(param, j, rhs.g0fun);
        end
    end
    
    function rval = feval(rhs, z)
        sz = size(z);
        z = z(:);
        rval = nan(numel(z), 1);
        [d, q, m] = domainDataB(rhs.domain);
%         alpha = rhs.parameter;
        
        for j = 0:m
            onCj = abs(q(j+1) - abs(z - d(j+1))) < eps(2);
            if ~any(onCj)
                continue
            end
            zj = z(onCj);
            
            if j == 0
                rval(onCj) = 2*pi*real(rhs.g0fun.hat(zj)) ...
                    + unwrap(angle(rhs.g0fun.singCorrFact(zj)));
                continue
            end
            rval(onCj) = 2*pi*real(rhs.gjFuns{j}.hat(zj));
%             thja = d(j+1) + q(j+1)^2/conj(alpha - d(j+1));
%             rval(onCj) = 2*pi*real(rhs.gjFuns{j}.hat(zj)) ...
%                 + unwrap(angle( (zj - d(j+1)).*(alpha - d(j+1)) ...
%                 ./(zj - alpha)./(zj - thja) ));
        end
        rval = reshape(rval, sz);
    end
end

end
