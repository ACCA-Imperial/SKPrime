%% Prime function usage example.

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
% test line

clear

%% Create a bounded circular domain and pick a parameter point.

dv = [0.0524476+0.365385i; -0.27972-0.127622i; 0.482517-0.281469i];
qv = [0.151967; 0.179551; 0.209557];
m = numel(dv);

alpha = -0.5i;


%% Compute the prime function (and its square).
% Note X is just the square of w, but using |X| skips some root branch
% checking done when evaluating the prime function.

w = skprime(alpha, dv, qv);
X = @(z) w.X(z);


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
% Note that calling the |invParam| method (from the |skprime| class) is
% equivalent to calling
%
%    w2 = skprime(1/conj(alpha), w);
%
% where using the previous instance of |skprime| accelerates the
% construction of the next (e.g., the first-kind integral functions do not
% need to be recomputed). Both methods are faster than calling
%
%    w2 = skprime(1/conj(alpha), dv, qv);

wi = invParam(w);
W = @(z) log(w(z)./(abs(alpha)*wi(z)))/(2i*pi);
Wp = W(zp);


%%
% Slightly more complicated potential; sum of 2 logarithmic singularities.

a2 = 0.50789+0.29737i;
w2 = skprime(a2, w);
w2i = invParam(w2);
W2 = @(z) W(z) + log(w2(z)./(abs(a2)*w2i(z)))/(2i*pi);


%%
% Plot some equipotential lines.

% Grid points in domain.
[X, Y] = meshgrid(linspace(-1, 1, 200));
zg = complex(X, Y);
zg(abs(zg) >= 1-eps(2)) = nan;
for j = 1:m
    zg(abs(zg - dv(j)) < qv(j)+eps(2)) = nan;
end

% Boundary points.
zb = exp(2i*pi*(0:200)'/200);
zb = bsxfun(@plus, [0; dv].', bsxfun(@times, [1; qv]', zb));

clf
contour(real(zg), imag(zg), imag(W2(zg)), 20, 'color', lines(1))
hold on
plot(zb, 'k-', 'linewidth', 1.5)
hold off
set(gca, 'dataaspectratio', [1, 1, 1])
axis off
