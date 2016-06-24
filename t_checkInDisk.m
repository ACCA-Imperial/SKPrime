%% Check "in disk" version.
clear

j = 1;
L = 7;


%%

dv = [
  -0.2517+0.3129i
   0.2307-0.4667i];
qv = [
  0.2377
  0.1557];
D = skpDomain(dv, qv);

% alpha = 0.53528 + 0.17493i;
prealpha = (-0.66122+0.28688i);
% alpha = dv(1) + qv(1);

theta = @D.theta;

alpha = theta(j, 1/conj(prealpha)); % <-------- the parameter of concern
% alpha = 1/conj(alpha);


%%

wja = skpindisk(alpha, D);


%%
% Check using product formula.

wref = skprod(dv, qv, L);

% figure(1), clf
% plotfd(@(z) wref(z, thja)./wref(z, 1/conj(alpha)), D, 'auto')
% figure(2), clf
% plotfd(rootHj, D, 'auto')

zp = 0.5;
err = wref(zp, alpha) - wja(zp);
disp(abs(err))
