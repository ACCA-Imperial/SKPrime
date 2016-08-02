classdef SingularityCorrectionGj < SKP.evaluable
%SingularityCorrectionGj encapsulates the singularity correction factor for
%the Greens function wrt Cj.
%
%See also greensCj.

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
    function sf = SingularityCorrectionGj(alpha, j, domain)
        if ~nargin
            return
        end
        
        % FIXME: validate input.
        %   alpha should be a skpParameter object.
        %   domain should be a skpDomain object.
        
        d = domain.centers;
        theta = @domain.theta;
        ak = isclose(domain, alpha);
        for k = ak(ak ~= j)
            if k > 0
                cfn = @(z) z - d(k);
                cfd = @(z) z - theta(k, 1/conj(alpha));
            else
                cfn = @(z) 1;
                cfd = @(z) z - 1/conj(alpha);
            end
            sf.correctionNumerators{end+1} = cfn;
            sf.correctionDenominators{end+1} = cfd;
        end
    end
    
    function dsf = diff(sf, n)
        %order n variable derivative.
        
        if isempty(sf.correctionNumerators)
            dsf = @(z) complex(zeros(size(z)));
            return
        end
        
        if nargin < 2
            n = 1;
        end
        % FIXME: validate input n.
        
        dmult = (-1)^(n-1)*factorial(n-1);
        
        function dv = deval(z)
            dv = 0;
            cfn = sf.correctionNumerators;
            cfd = sf.correctionDenominators;
            for i = 1:numel(cfn)
                dv = dv + dmult*(1./cfn{i}(z).^n - 1./cfd{i}(z).^n);
            end
        end
        
        dsf = @deval;
    end
    
    function v = feval(sf, z)
        %provides function evaluation behaviour.
        
        v = 1;
        cfn = sf.correctionNumerators;
        cfd = sf.correctionDenominators;
        for i = 1:numel(cfn)
            v = v.*cfn{i}(z)./cfd{i}(z);
        end
    end
end

end
