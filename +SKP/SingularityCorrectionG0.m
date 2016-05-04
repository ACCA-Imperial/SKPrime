classdef SingularityCorrectionG0 < evaluable
%SingularityCorrectionG0 encapsulates the singularity correction factor for
%the Greens function wrt C0.
%
%See also greensC0.

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

properties(SetAccess=protected)
    correctionNumerators = {}
    correctionDenominators = {}
end

methods
    function sf = SingularityCorrectionG0(alpha, domain)
        if ~nargin
            return
        end
        
        % FIXME: validate input.
        %   alpha should be a skpParameter object.
        %   domain should be a skpDomain object.
        
        [d, q] = domainData(domain);
        aj = isclose(domain, alpha);
        for j = aj
            if alpha == 0
                thj0 = d(j);
                thjinf = d(j) - q(j)^2/conj(d(j));
                cfn = @(z) (z - thj0);
                cfd = @(z) (z - thjinf);
            else
                thj = @(z) d(j) + q(j)^2*z./(1 - conj(d(j))*z);
                cfn = @(z) (z - thj(alpha));
                cfd = @(z) (z - thj(inv(alpha)));
            end
            sf.correctionNumerators{end+1} = cfn;
            sf.correctionDenominators{end+1} = cfd;
        end
    end
    
    function dsf = diff(sf, n)
        %Order n variable derivative.
        
        if nargin < 2
            n = 1;
        end
        % FIXME: validate input n.
        
        dmult = (-1)^(n-1)*factorial(n-1);
        
        function dv = deval(z)
            dv = complex(zeros(size(z)));
            
            cfn = sf.correctionNumerators;
            cfd = sf.correctionDenominators;
            for i = 1:numel(cfn)
                dv = dv + dmult*(1./cfn{i}(z) - 1./cfd{i}(z));
            end
        end
        
        dsf = @deval;
    end
    
    function v = feval(sf, z)
        %Provide function evaluation behaviour.
        
        v = complex(ones(size(z)));
        cfn = sf.correctionNumerators;
        cfd = sf.correctionDenominators;
        for i = 1:numel(cfn)
            v = v.*cfn{i}(z)./cfd{i}(z);
        end
    end
end

end
