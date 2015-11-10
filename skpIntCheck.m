function rerr = skpIntCheck(f, z0, acal, D, N, Nf)
%skpIntCheck checks the skprime function via integral formula.
%
% rerr = skpIntCheck(f, z0, acal, D, N, Nf)
%   f is the function handle for comparison.
%   z0 is a vector of points to check.
%   acal is the point to use for "calibration" of the integral formula,
%     (point used to compute the constant) -- the point just must be in
%     the domain.
%   D is a circle domain (skpDomain object).
%   N [optional] is the number of spectral points.
%   Nf [optional] is the number of Fourier points to use.
% The function returns an array, size(z0), of relative error measured
% at z0.
%
% The function checks the accuracy of the computed Schottky-Klein
% prime function by using a known comparision function and the
% formula from [D. Crowdy, "The Schwarz problem in multiply
% connected domains and the Schottky-Klein prime function", Complex
% Variables and Elliptic Equations, 53(3) 2008, 221-236],
% 
%                        /
%                   1    |
%   g(a) = 1i*C + -----  |  real(f(z)) * [d log(w(z,a)) ...
%                 2i*pi  |                  + conj(d log(w(z,1/conj(a)))]
%                        / 
%                      boundary(D)
%
% where w(z,a) is the Shottky-Klein prime function and C is a constant
% which must be computed. The DFT with an even number of Fourier points
% is used to compute the derivative of the prime function on the
% boundary, while the integral is computed via the trapezoid rule. Then
%
%  rerr = abs(f(a) - g(a))/abs(f(a))
%
% is the relative error.
%
% See also: skprime, skpDomain

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

if nargin < 6
    Nf = 256;
elseif mod(Nf, 2)
    error('Number of Fourier points must be a multiple of 2.')
end

if nargin < 5
    N = 32;
end

accel = false;
if iscell(D) && all(cellfun(@(c) isa(c, 'vjFirstKind'), D(:)))
    accel = true;
    vjFuns = D;
    D = vjFuns{1}.domain;
    Nv = vjFuns{1}.truncation;
    N = max(N, Nv);
end

% Boundary points
[d, q, m] = domainDataB(D);
eit = exp(2i*pi/Nf*(0:Nf-1)');
zj = bsxfun(@plus, d.', bsxfun(@times, q', eit));

% Real part of f along with boundary orientation.
sig = @(j) 1 - (j > 1)*2;
phi = bsxfun(@times, sig(1:m+1), real(f(zj)));

% Derivative on boundary via DFT.
function dfdt = dftDeriv(fz)
    dmult = [0:Nf/2-1, 0, -Nf/2+1:-1]';
    dfdt = complex(nan(size(fz)));
    for j = 1:m+1
        dfdt(:,j) = ifft(1i*dmult.*fft(fz(:,j)));
    end
end

% Integral formula.
function I = skIntegral(w, wc)
    wz = w(zj);
    dlogw = dftDeriv(wz)./wz;
    wcz = wc(zj);
    dlogwc = dftDeriv(wcz)./wcz;
    
    I = sum(phi(:).*(dlogw(:) + conj(dlogwc(:))))/(1i*Nf);
end


% Calibration (compute constant)
if accel
    w = skprime(acal, vjFuns, N);
else
    w = skprime(acal, D, N);
end
wc = skprime(1/conj(acal), w);

I = skIntegral(w, wc);
C = real(-1i*(f(acal) - I));


% Measure error at each point.
fhat = complex(nan(size(z0)));
for k = 1:numel(z0)
    w = skprime(z0(k), w, N);
    wc = skprime(1/conj(z0(k)), w, N);
    
    I = skIntegral(w, wc);
    fhat(k) = I + 1i*C;
end
rerr = abs(f(z0) - fhat)./abs(f(z0));

end
