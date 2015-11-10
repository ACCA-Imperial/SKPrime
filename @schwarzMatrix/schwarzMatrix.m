classdef schwarzMatrix < handle
%schwarzMatrix provides shared memory for domain matrix.

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

properties
    matrix
end

methods
    function M = schwarzMatrix(mat)
        if ~nargin
            return
        end
        
        M.matrix = mat;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function C = plus(A, B)
        C = A.matrix + B;
    end
    
    function C = uplus(A)
        C = +A.matrix;
    end
    
    function C = minus(A, B)
        C = A.matrix - B;
    end
    
    function C = uminus(A)
        C = -A.matrix;
    end
    
    function C = times(A, B)
        C = A.matrix.*B;
    end
    
    function C = rdivide(A, B)
        C = A.matrix./B;
    end
    
    function C = ldivide(A, B)
        C = A.matrix.\B;
    end
    
    function C = power(A, B)
        C = A.matrix.^B;
    end
    
    function C = mtimes(A, B)
        C = A.matrix*B;
    end
    
    function C = mrdivide(A, B)
        C = A.matrix/B;
    end
    
    function C = mldivide(A, B)
        C = A.matrix\B;
    end
    
    function C = mpower(A, B)
        C = A.matrix^B;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function A = subsasgn(A, S, B)
        A.matrix = builtin('subsasgn', A.matrix, S, B);
    end
    
    function B = subsref(A, S)
        if S.type == '.' && strcmp(S.subs, 'matrix')
            B = A.matrix;
            return
        end
        B = builtin('subsref', A.matrix, S);
    end
end

end
