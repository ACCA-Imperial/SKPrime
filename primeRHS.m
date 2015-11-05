classdef primeRHS < skpObject
%primeRHS computes RHS for Schwarz problem for finding Prime function.

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
    
    vjFuns
    g0fun
    prfact
end

properties(Access=private)
    bvfun
end

methods
    function rhs = primeRHS(param, vjfuns)
        if ~nargin
            return
        end
        
        if ~isa(param, 'skpParameter')
            error('SKPrime:invalidArgument', ...
                'Parameter must be a "skpParameter" object.')
        end
        rhs.parameter = param;
        rhs.domain = vjfuns{1}.domain;
        rhs.vjFuns = vjfuns;

        if ~isempty(param.ison)
            rhs.bvfun = @rhs.onBdryFun;
            return
        end

        rhs.g0fun = G0Cauchy(param, vjfuns{1});

        % Is parameter near outer boundary?
        [d, q] = domainData(rhs.domain);
        thf = @(z) 1;
        for j = isclose(rhs.domain, param)'
            thj = @(z) d(j) + q(j)^2*z./(1 - conj(d(j))*z);
            thf = @(z) thf(z).*(z - thj(z)).^2 ...
                ./(z - thj(param))./(z - thj(inv(param)));
        end
        rhs.prfact = thf;

        rhs.bvfun = @rhs.inDomainFun;
    end
    
    function rval = feval(rhs, z)
        sz = size(z);
        z = z(:);
        rval = nan(prod(sz), 1);
        [d, q, m] = domainDataB(rhs.domain);
        
        for j = 0:m
            onCj = abs(q(j+1) - abs(z - d(j+1))) < eps(2);
            if ~any(onCj)
                continue
            end
            
            rval(onCj) = rhs.bvfun(j, z(onCj));
        end
        rval = reshape(real(rval), sz);
    end
end

methods(Access=protected)
    function val = inDomainFun(rhs, j, zj)
        val = 2*pi*real(rhs.g0fun.hat(zj));
        if j > 0
            vj = rhs.vjFuns{j};
            val = val + ...
                2*pi*real(vj(rhs.parameter) - vj.hat(zj)) + Aj(rhs, j, zj);
        else
            val = val + A0(rhs, zj);
        end
    end
    
    function val = onBdryFun(rhs, j, zj)
        alpha = rhs.parameter;
        ja = alpha.ison;
        
        if j == ja
            val = complex(zeros(size(zj)));
            return
        end
        if j == 0
            va = rhs.vjFuns{ja};
            val = 2*pi*real(va.hat(zj) - va(alpha)) + A0ja(rhs, ja, zj);
        else
            vj = rhs.vjFuns{j};
            if ja == 0
                val = 2*pi*real(vj(alpha) - vj.hat(zj)) ...
                    + Aj0(rhs, j, zj);
            else
                va = rhs.vjFuns{ja};
                val = 2*pi*real(va.hat(zj) - vj.hat(zj) ...
                    + vj(alpha) - va(alpha)) + Ajja(rhs, j, ja, zj);
            end
        end
    end
    
    function val = A0(rhs, zj)
        alpha = rhs.parameter;
        thf = rhs.prfact;

        if abs(alpha) == 0 || isinf(alpha)
            val = thf(zj);
        else
            val = unwrap(angle( alpha*zj.*thf(zj) ...
                ./(zj - alpha)./(zj - inv(alpha)) ));
        end
    end
    
    function val = Aj(rhs, j, zj)
        alpha = rhs.parameter;
        [d, q, ~, di] = domainData(rhs.domain);
        thf = rhs.prfact;
        
        if abs(d(j)) > q(j)
            if 0 < abs(alpha) && ~isinf(alpha)
                val = unwrap(angle(...
                    alpha.*(zj - di(j)).*thf(zj) ...
                    ./(zj - alpha)./(zj - 1/conj(alpha)) ));
            else
                val = unwrap(angle( (zj - di(j))./zj).*thf(zj) );
            end
        else
            val = unwrap(angle(...
                alpha.*thf(zj)./(zj - alpha)./(zj - 1/conj(alpha)) ));
        end
    end
    
    function val = A0ja(rhs, ja, zj)
        alpha = rhs.parameter;
        [d, q, ~, di] = domainDataB(rhs.domain);
        dRatio = @(z,j) (z - d(j+1))./(z - di(j));
        
        if abs(d(ja+1)) > q(ja+1)
            val = unwrap(angle(...
                zj*(alpha - d(ja+1)).*dRatio(zj, ja)./(zj - alpha).^2));
        elseif abs(abs(d(ja+1)) - q(ja+1)) < eps(2) ...
                && (abs(alpha) == 0 || isinf(alpha))
            val = unwrap(angle(-d(ja+1)*(zj - d(ja+1))./zj));
        else
            val = unwrap(angle(...
                zj*(alpha - d(ja+1)).*(zj - d(ja+1))./(zj - alpha).^2));
        end
    end
    
    function val = Aj0(rhs, j, zj)
        alpha = rhs.parameter;
        [d, q, ~, di] = domainDataB(rhs.domain);
        
        if abs(d(j+1)) > q(j+1)
            val = unwrap(angle(alpha*(zj - di(j))./(zj - alpha).^2));
        else
            val = unwrap(angle(alpha./(zj - alpha).^2));
        end
    end
    
    function val = Ajja(rhs, j, ja, zj)
        alpha = rhs.parameter;
        [d, q, ~, di] = domainDataB(rhs.domain);
        dRatio = @(z,j) (z - d(j+1))./(z - di(j));
        
        daj = d(ja+1);
        dk = d(j+1);
        if abs(daj) > q(ja+1) && abs(dk) > q(j+1)
            val = unwrap(angle(...
                (alpha - daj)*(zj - di(j)).*dRatio(zj, ja) ...
                ./(zj - alpha).^2 ));
        elseif abs(daj) > q(ja+1)
            val = unwrap(angle(...
                (alpha - daj).*dRatio(zj, ja)./(zj - alpha).^2));
        else
            if 0 < abs(alpha) && ~isinf(alpha)
                val = unwrap(angle(...
                    (alpha - daj).*(zj - di(j)).*(zj - daj)./(zj - alpha).^2));
            else
                val = unwrap(angle(-daj.*(zj - di(j)).*(zj - daj)./zj.^2));
            end
        end
    end    
end

methods(Hidden)
    function rhs = flipConstState(rhs)
        rhs.minusConstant = ~rhs.minusConstant;
    end
end

end % primeRHS
