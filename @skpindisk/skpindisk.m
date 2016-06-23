classdef skpindisk < skprime
%skpindisk is the SKPrime function for the parameter in a 1st level
%refelction of the domain.

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
    auxParameter
    indisk
    
    hejhal
    rootHejhal
end

methods
    function skp = skpindisk(alpha, dv, qv)
        % What I really want to do:
        %  1) Process arguments.
        %  2) Use domain and parameter to determine disk and auxiliary
        %  parameter.
        %  3) Use auxiliary parameter to call superclass constructor.
        
        args = {};
        if nargin == 2
            args = {alpha, dv};
        elseif nargin == 3
            args = {alpha, dv, qv};
        end
        
        skp = skp@skprime(args{:});
        if ~nargin
            return
        end
        
        D = skp.domain;
        alpha = skp.parameter;
        ja = alpha.indisk;
        if alpha.state == paramState.innerDisk
            beta = D.theta(ja, 1/conj(alpha));
        elseif alpha.state == paramState.outerDisk
            beta = D.theta(ja, alpha);
        else
            error('SKPrime:invalidArgument', ...
                'The parameter is not in an inner or outer disk.')
        end
        
        skp.auxParameter = 1/conj(beta);
    end
end

end






















