classdef(Abstract) evaluable
%evaluable provides the function evaluation protocol.
%
%If C is an instance of a class adopting the evaluable protocol, then the
%instance handles the syntax
%
%  v = C(z)
%
%Classes adopting this protocol must provide a method with the signature
%
%  v = feval(obj, z)
%
%which is called from the evaulable subsref.

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

methods
    v = feval(obj, z)
    
    function out = subsref(obj, S)
        %Provides w = f(z) syntax for object.
        
        if numel(S) == 1 && strcmp(S.type, '()')
            out = feval(obj, S.subs{:});
        else
            out = builtin('subsref', obj, S);
        end
    end
end

end
