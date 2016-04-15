classdef genericPlusSingular < bvpFun
%genericPlusSingular represents a generic BVP.
%
%For a single-valued imaginary function on a bounded, multiply connected
%domain D, we assume we know the imaginary part of the function up to an
%unknown constant. The modified Schwarz problem is given on the boundary of
%the domain by
%
%   f(z) = phi(z) + 1i*k + 1i*known(z)
%
%where the unknown constant k is in general a different value on each
%boundary.
%
%  g = genericPlusSingular(singPart, known, D)
%    Solves the modified Schwarz problem using known which is the function
%    handle to the imaginary part of the function on the boundary. The
%    final constructed function is given by
%       g(z) = f(z) + singPart(z)
%    where singPart is the handle to the singular part of the final
%    function.
%
%See also bvpFun.

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
    knownImaginary
    singularPart
end

methods
    function f = genericPlusSingular(singPart, known, D)
        if ~nargin
            args = {};
        else
            args = {D};
        end
        
        f = f@bvpFun(args{:});
        if ~nargin
            return
        end
        
        f.phiFun = solve(f.phiFun, known);
        f.boundaryFunction = @(z) f.phiFun(z) + 1i*known(z);
        f.continuedFunction = SKP.bmcCauchy(...
            f.boundaryFunction, f.domain, 2*f.truncation);
        
        f.knownImaginary = known;
        f.singularPart = singPart;
    end
    
    function v = feval(f, z)
        %provides function evaluation.
        
        v = f.hat(z) + f.singularPart(z);
    end
    
    function v = hat(f, z)
        %evaluates the "analytic" part of the function.
        
        v = bvpEval(f, z);
    end
end

end
