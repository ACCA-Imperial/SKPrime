%% Integral prime check.
% Quick check of prime function accuracy. See help text in skpIntCheck for
% more information.
clear

dv = [0.5; -0.1+0.35i; -0.4i];
qv = 0.3*ones(size(dv));
D = skpDomain(dv, qv);

% Randomly chosen point in the domain to calibrate check integral.
acal = -0.3-0.1i;

% Test function of poles in holes.
f = @(z) reshape(sum(1./bsxfun(@minus, z(:), dv.'), 2), size(z));

% Some random points to check.
z0 = [0.4+0.5i; -0.5-0.5i; +0.5-0.5i; -0.5+0.1i; 0.8i];

tic
rerr = skpIntCheck(f, z0, acal, D);
toc
disp(rerr)
