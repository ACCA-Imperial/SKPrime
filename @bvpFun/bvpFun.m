classdef bvpFun < skpObject
%bvpFun is the base class for the BVP functions.
%
% obj = bvpFun(D, N, phi)
%   Base class for SKPrime BVP based functions. Not a truly functional
%   class.
%
% Provides properties:
%   phiFun
%   trunctation
%
% See also skpObject.

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
    phiFun                      % Unknown function object
    truncation = 64             % Series truncation setting
end

methods
    function bvp = bvpFun(D, N, phi)
        if ~nargin
            sarg = {};
        elseif isa(D, 'bvpFun')
            phi = D.phiFun;
            if nargin < 2
                N = D.truncation;
            end
            sarg = {D.domain};
        else
            sarg = {D};
            if nargin < 2
                N = [];
            end
            if nargin < 3
                phi = [];
            end
        end
        bvp = bvp@skpObject(sarg{:});
        if ~nargin
            return
        end
        
        if ~isempty(N)
            if ~isnumeric(N) || N <= 0 || N ~= floor(N)
                error('SKPrime:invalidArgument', ...
                    'N must be a non-zero, positive integer.')
            end
            bvp.truncation = N;
        end
        
        if isempty(phi)
            bvp.phiFun = schwarz(bvp.domain, bvp.truncation);
        else
            if ~isa(phi, 'schwarz')
                error('SKPrime:invalidArgument', ...
                    '"phi" must be a "schwarz" object.')
            end
            bvp.phiFun = phi;
        end
    end
end

end
