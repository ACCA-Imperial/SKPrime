%% Hejhal identity temp code.
clear

j = 1;


%%

dv = [
  -0.2517+0.3129i
   0.2307-0.4667i];
qv = [
  0.2377
  0.1557];
D = skpDomain(dv, qv);

% alpha = 0.53528 + 0.17493i;
alpha = -0.66122+0.28688i;
% alpha = dv(1) + qv(1);

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);

thja = theta(j, 1/conj(alpha)); % <-------- the parameter of concern


%%
% Step 1: Identify which circle the parameter is in.

try
    thja = skpParameter(thja, D);
catch err
    fprintf('%s\n', err.message)
end

if isa(thja, 'skpParameter') && thja.state == paramState.onInnerBdry
    fprintf('On inner boundary %d. Hejhal not needed.\n', thja.ison)
    return
end

ja = find(abs(thja - dv) <= qv);
fprintf('Parameter is inside circle %d.\n', ja);


%%
% Step 2: Define the "pair parameter" beta in the domain, and find the prime
% function for this parameter.

beta = 1/conj(theta(ja, 1/conj(thja)));

if ~(isin(D, beta) || isin(D, 1/conj(beta)))
    fprintf('Parameter is more than one Schottky group level deep.\n')
    return
end
if abs(1/conj(alpha) - beta) > eps(20)
    fprintf('Beta has unexpected value.\n')
    return
end

wb = skprime(beta, dv, qv);


%%
% Step 3: Use this to define the prime function for alpha via the Hejhal
% transformation.

vj = wb.vjFuns{ja};
taujj = vj.taujj;
rootHj = @(z) -exp(-2i*pi*(vj(beta) - vj(z) + taujj/2)) ...
    *qv(ja)/(1 - conj(dv(j))*beta);

wja = @(z) wb(z).*rootHj(z);

% clf
% plotfd(wja, D, 'auto')


%%
% Check using product formula.

wref = skprod(dv, qv, 6);

% figure(1), clf
% plotfd(@(z) wref(z, thja)./wref(z, 1/conj(alpha)), D, 'auto')
% figure(2), clf
% plotfd(rootHj, D, 'auto')

zp = 0.5;
err = wref(zp, theta(j, 1/conj(alpha))) - wja(zp);
disp(err)


%%
% Check using Gj?

gjref = greensCj(alpha, j, wb);
egjref = @(z) exp(2i*pi*gjref(z));

wa = skprime(alpha, wb);
egjan = @(z) wa(z)./wja(z)*qv(j)/abs(alpha - dv(j));
gjan = @(z) log(egjan(z))/2i/pi;

% figure(1), clf
% plotfd(gjref, D, 'unit')
% figure(2), clf
% plotfd(gjan, D, 'unit')

zp = 0.5;
err = gjref(zp) - gjan(zp);
disp(err)






























