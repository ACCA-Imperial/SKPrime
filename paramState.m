classdef paramState < int8
%paramState describes the state of the alpha parameter.

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

enumeration
    isZero(-3)
    innerFD(-2)
    onInnerBdry(-1)
    isUnit(0)
    onOuterBdry(1)
    outerFD(2)
    atInf(3)
end

methods
    function state = inv(state)
        m = enumeration(state);
        state = m(-state + (numel(m) + 1)/2);
    end
end

end
