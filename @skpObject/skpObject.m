classdef skpObject
%skpObject abstract base utility class for SKPrime
%
% obj = skpObject(domain)
%   Abstract base utility class for the SKPrime kit. Must be subclassed.
%
% Abstract method:
%   feval -- define to provide function evaluation behaviour

% Everett Kropf, 2015
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
    domain                      % skpDomain object
end

methods(Abstract)
    v = feval(skp, z)
end

methods
    function skp = skpObject(D)
        if ~nargin
            return
        end
        
        if isa(D, 'circleRegion')
            D = skpDomain(D);
        end
        
        if ~isa(D, 'skpDomain')
            error('SKPrime:invalidArgument', 'Expected a "skpDomain" object.');
        end
        skp.domain = D;
    end
    
    function onB = onBoundary(skp, z)
        %checks for points on boundary
        %
        % onB = onBoundary(skp, z)
        %   Return boolean array of size(z) where true if a point z is on
        %   the boundary, false otherwise.
        
        onB = false(size(z));
        [d, q, m] = domainDataB(skp.domain);
        for j = 0:m
            onB(abs(q(j+1) - abs(z - d(j+1))) < eps(2)) = true;
        end
    end
    
    function out = subsref(skp, S)
        %provides function syntax
        %
        % Provides skp(z) evaluation syntax, passes all unrecognized
        % subsref to builtin.
        
        if numel(S) == 1 && strcmp(S.type, '()')
            out = feval(skp, S.subs{:});
        else
            out = builtin('subsref', skp, S);
        end
    end
end

end
