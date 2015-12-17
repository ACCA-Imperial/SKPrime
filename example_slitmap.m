%% Circles to horizontal slits via prime function.
% Map from a bounded circle domain to a horizontal slit domain.

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

dv = [
    -0.2517+0.3129i
    0.2307-0.4667i ];
qv = [
    0.2377
    0.1557 ];
m = numel(dv);

alpha = -0.21789-0.22421i;
beta = 0.52105+0.20526i;

nb = 200;
zb = boundaryPts(skpDomain(dv, qv), nb);
zb = zb(:,1:m+1);
zb(end+1,:) = zb(1,:);


%%
% Circular slit map.

w1 = skprime(alpha, dv, qv);
w1i = invParam(w1);
w2 = skprime(beta, w1);
w2i = invParam(w2);

P = @(z) w1(z).*w2i(z)./w2(z)./w1i(z);

Pzb = P(zb);


%%
% Log map for vertical slits, then rotation.

Hzb = exp(-1i*pi/2)*log(Pzb);


%%

aspecteq = @() set(gca, 'dataaspectratio', [1, 1, 1]);

figure(1), clf


figure(3), clf
plot(Hzb)
aspecteq()
axis off
