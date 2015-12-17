%% Integral prime check.
% Quick check of prime function accuracy. See help text in skpIntCheck for
% more information.

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

clear

dv = [0.5; -0.1+0.35i; -0.4i];
qv = 0.3*ones(size(dv));

% Randomly chosen point in the domain to calibrate check integral.
acal = -0.3-0.1i;

% Test function of poles in holes.
f = @(z) reshape(sum(1./bsxfun(@minus, z(:), dv.'), 2), size(z));

% Some random points to check.
z0 = [0.4+0.5i; -0.5-0.5i; +0.5-0.5i; -0.5+0.1i; 0.8i];

rerr = skpIntCheck(f, z0, acal, skpDomain(dv, qv));


%%

fprintf('\nRelative error at the test points:\n\n')
fprintf('    z0       rel.err.\n')
fprintf(' ---------  ----------\n')
for i = 1:numel(z0)
    fprintf(' %9s  %.5g\n', num2str(z0(i)), rerr(i))
end
fprintf('\n')
