%% Hejhal identity temp code.
clear

j = 1;
L = 6;


%%

dv = [
  -0.2517+0.3129i
   0.2307-0.4667i];
qv = [
  0.2377
  0.1557];
D = skpDomain(dv, qv);

% alpha = 0.53528 + 0.17493i;
prealpha = 1/conj(-0.66122+0.28688i);
% alpha = dv(1) + qv(1);

theta = @D.theta;

alpha = theta(j, 1/conj(prealpha)); % <-------- the parameter of concern
% alpha = 1/conj(alpha);


%%
% Step 1: Identify which circle the parameter is in and define the "pair
% parameter" in the FD.

try
    alpha = skpParameter(alpha, D);
catch err
    fprintf('%s\n', err.message)
end

if isa(alpha, 'skpParameter') && alpha.state == paramState.onInnerBdry
    fprintf('On inner boundary %d. Hejhal not needed.\n', alpha.ison)
    return
end

if abs(alpha) < 1
    ja = find(abs(alpha - dv) < qv + eps(2), 1);
    if ~isempty(ja)
        beta = theta(ja, 1/conj(alpha));
        bcond = (isin(D, beta) || isin(D, 1/conj(beta)));
    end
    if isempty(ja) || ~bcond
        fprintf('Parameter is more than one Schottky group level deep.\n')
        return
    end
    fprintf('Parameter is in 1st FD copy inside circle %d.\n', ja);
else
    ja = find(abs(1/conj(alpha) - dv) < qv + eps(2), 1);
    if ~isempty(ja)
        beta = theta(ja, alpha);
        bcond = isin(D, beta) || isin(D, 1/conj(beta));
    end
    if isempty(ja) || ~bcond
        fprintf('Parameter is more than one Schottky group level deep.\n')
        return
    end
    fprintf('Parameter is in 1st FD copy inside circle %d''.\n', ja);
end


%%
% Step 3: Use the secondary parameter to define the prime function for the
% main parameter.

beta = 1/conj(beta);
wb = skprime(beta, dv ,qv);
vj = wb.vjFuns{ja};
taujj = vj.taujj;
rootHj = @(z) -exp(-2i*pi*(vj(beta) - vj(z) + taujj/2)) ...
    *qv(ja)/(1 - conj(dv(j))*beta);

if abs(alpha) < 1
    wja = @(z) wb(z).*rootHj(z);
else
    wja = @(z) -(z/conj(theta(ja, beta))) ...
        .*conj(wb(1./conj(z)).*rootHj(1./conj(z)));
end

% clf
% plotfd(wja, D, 'auto')


%%
% Check using product formula.

wref = skprod(dv, qv, L);

% figure(1), clf
% plotfd(@(z) wref(z, thja)./wref(z, 1/conj(alpha)), D, 'auto')
% figure(2), clf
% plotfd(rootHj, D, 'auto')

zp = 0.5;
err = wref(zp, alpha) - wja(zp);
disp(err)


%%
% Check using Gj?

% gjref = greensCj(alpha, j, wb);
% egjref = @(z) exp(2i*pi*gjref(z));
% 
% wa = skprime(alpha, wb);
% egjan = @(z) wa(z)./wja(z)*qv(j)/abs(alpha - dv(j));
% gjan = @(z) log(egjan(z))/2i/pi;
% 
% % figure(1), clf
% % plotfd(gjref, D, 'unit')
% % figure(2), clf
% % plotfd(gjan, D, 'unit')
% 
% zp = 0.5;
% err = gjref(zp) - gjan(zp);
% disp(err)






























