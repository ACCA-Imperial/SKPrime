%% Circles to horizontal slits via prime function.
% Map from a bounded circle domain to a horizontal slit domain. The
% circular slit map given below is derived in terms of the Schotty-Klein
% prime function in
%
%    D. Crowdy and J. Marshall. "Conformal mappings between canonical
%    multiply connected domains." CMFT, 6(1), 2006, pp. 59--76.

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


%%
% Some domain setup.

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

% Points on the boundary.
zb = exp(2i*pi/200*(0:200)');
zb = bsxfun(@plus, [0; dv].', bsxfun(@times, [1; qv]', zb));


%%
% Some setup for later plots.

% Make some grid points in the circle domain.
np = 80;
nlines = 15;
zh = bsxfun(@plus, 1i*linspace(-1, 1, nlines), linspace(-1, 1, np)');
zv = bsxfun(@plus, linspace(-1, 1, nlines), 1i*linspace(-1, 1, np)');
zh(abs(zh) >= 1-eps(2)) = nan;
zv(abs(zv) >= 1-eps(2)) = nan;
for j = 1:m
    zh(abs(zh - dv(j)) < qv(j)+eps(2)) = nan;
    zv(abs(zv - dv(j)) < qv(j)+eps(2)) = nan;
end

% Plot helpers.
aspecteq = @() set(gca, 'dataaspectratio', [1, 1, 1]);
zcolor = lines(2);
plotv = @(z) plot(z, '.', 'color', zcolor(1,:));
ploth = @(z) plot(z, '.', 'color', zcolor(2,:));
plotb = @(z) plot(z, 'k-', 'linewidth', 1.5);


%%
% Plot the circle domain for reference.

figure(1), clf
hold on
plotv(zv)
ploth(zh)
plotb(zb)
hold off
aspecteq()
axis off


%%
% Construct a circular slitmap. This is a map from the bounded circle
% domain to the domain with concentric circular arc slits with respect to
% the origin. Let |w| represent the SK-prime function. Then the map takes
% the form
%
%           w(z, alpha) * w(z, 1/conj(beta))
%    P(z) = --------------------------------
%           w(z, beta) * w(z, 1/conj(alpha))

w1 = skprime(alpha, dv, qv);
w1i = invParam(w1);
w2 = skprime(beta, w1);
w2i = invParam(w2);

P = @(z) w1(z).*w2i(z)./w2(z)./w1i(z);


%%
% A plot of the image of the circle domain under P(z).

figure(2), clf
hold on
plotv(P(zv))
ploth(P(zh))
plotb(P(zb))
hold off
aspecteq()
axis(0.8*[-1, 1, -1, 1])
axis off


%%
% Log map for vertical slits, then rotation.

H = @(z) exp(-1i*pi/2)*log(P(z));
Hzb = H(zb);


%%
% A plot of the circle domain under H(z).

figure(3), clf
hold on
plotv(H(zv))
ploth(H(zh))
plotb(Hzb)
hold off
aspecteq()
axis off
