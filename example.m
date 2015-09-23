%% Prime function usage example.
% E. Kropf, 2015

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


%% Create a bounded circular domain and pick a parameter point.

dv = [0.0524476+0.365385i; -0.27972-0.127622i; 0.482517-0.281469i];
qv = [0.151967; 0.179551; 0.209557];
m = numel(dv);
D = skpDomain(dv, qv);

alpha = -0.5i;


%% Compute the prime function (and its square).
% Note X is just the square of w, but using |Xeval| skips some root branch
% checking done when evaluating the prime function.

w = skprime(alpha, D);
X = @(z) Xeval(w, z);


%% Evalutate the functions at various points.

np = 10;
zp = complex(nan(np, m));
for j = 1:m
    zp(:,j) = dv(j) + 1.2*qv(j)*exp(2i*pi*(0:np-1)'/np);
end
zp = [zp, 1./conj(zp)];

wp = w(zp);
Xp = X(zp);


%% Make a complex potential.
% Note using the first prime function to construct the second reduces the
% amount of work done (the v_j functions do not need recomputed, and may be
% copied).

w2 = skprime(1./conj(alpha), w);
W = @(z) log(w(z)./(abs(alpha)*w2(z)))/(2i*pi);
Wp = W(zp);
