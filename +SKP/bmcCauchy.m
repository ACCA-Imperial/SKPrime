function f = bmcCauchy(h, D, N)
%bmcCauchy is the Cauchy interpolant for bounded MC circle domain.
%
% f = bmcCauchy(h, D, N)
%   h is boundary function
%   D is circle domain
%   N is number of points in collocation

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

if isa(D, 'circleRegion')
    D = skpDomain(D);
end
[d, q, m] = domainDataB(D);

% Collocation points and boundary data.
sig = @(j) 1 - (j > 0)*2;
eit = exp(2i*pi*(0:N-1)'/N);
tjk = zeros((m+1)*N, 1);
sjk = tjk;
hjk = tjk;
for j = 1:m+1
    k = (j-1)*N+(1:N);
    sjk(k) = sig(j-1)*q(j)*eit;
    tjk(k) = d(j) + sig(j-1)*sjk(k);
    hjk(k) = h(tjk(k)).*sjk(k);
end

nsrc = numel(tjk);
tjkv = [real(tjk), imag(tjk)].';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = feval(z)
    persistent hasfmm
    if isempty(hasfmm)
        if exist('zfmm2dpart', 'file') == 2
            hasfmm = true;
        else
            hasfmm = false;
        end
    end
    
    if hasfmm && numel(z) > 1000
        f = fevalFmm(z);
    else
        f = fevalMat(z);
    end
end

function f = fevalMat(z)
    f = complex(nan(size(z)));
    I = 1./bsxfun(@minus, tjk.', z(:));
    f(:) = (I*hjk)./(I*sjk);
end

function f = fevalFmm(z)
    f = complex(nan(size(z)));    
    notnan = ~isnan(z);
    ntarg = sum(notnan(:));
    z = [real(z(notnan)), imag(z(notnan))].';
    
    I = zfmm2dpart(5, nsrc, tjkv, -hjk, false, false, false, ...
        ntarg, z, true, false, false);
    J = zfmm2dpart(5, nsrc, tjkv, -sjk, false, false, false, ...
        ntarg, z, true, false, false);
    f(notnan) = I.pottarg./J.pottarg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = @feval;

end % bmcCauchy
