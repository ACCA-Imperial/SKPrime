function install(arg)
%install puts SKPrime in search path.
%
% Installation is defined simply to be the addition of the SKPrime
% directory to MATLAB's search path. Uninstall removes this directory from
% the search path.
%
% install
%   performs installation operations for SKPrime.
% install -u
%   performs uninstall operations for SKPrime.

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

if exist('./@schwarz/schwarz.m', 'file') ~= 2 || ...
        exist('./@skpfunction/skpfunction.m', 'file') ~= 2 || ...
        exist('./@skpObject/skpObject.m', 'file') ~= 2 || ...
        exist('./@paramState/paramState.m', 'file') ~= 2
    error('skprime:RuntimeError', ...
        'This utility must be run from the SKPrime directory.')
end
skpdir = pwd;

if nargin && strcmp(arg, '-u')
    uninstall(skpdir)
    return
end

fprintf('Adding directory "%s" to your MATLAB search path ... ', skpdir)
try
    addpath(skpdir)
    savepath
catch me
    fprintf('FAIL\n')
    rethrow(me)
end
fprintf('OK\n')

fprintf('SKPrime install complete. Try "example.m" now.\n')

end

function uninstall(skpdir)

fprintf('Removing directory "%s" from your MATLAB search path ... ', skpdir)
try
    rmpath(skpdir)
    savepath
catch me
    fprintf('FAIL\n')
    rethrow(me)
end
fprintf('OK\n')

fprintf('SKPrime uninstall complete. This directory may be safely deleted.\n')

end
