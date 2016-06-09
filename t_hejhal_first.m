%% Hejhal identity temp code.
clear

j = 2;


%%

dv = [
  -0.2517+0.3129i
   0.2307-0.4667i];
qv = [
  0.2377
  0.1557];
D = skpDomain(dv, qv);

alpha0 = 0.53528 + 0.17493i;
% alpha = dv(1) + qv(1);

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);

alpha = theta(j, 1/conj(alpha0)); % <-------- the parameter of concern


%%
% Step 1: Identify which circle the parameter is in.

try
    alpha = skpParameter(alpha, D);
catch err
    fprintf('%s\n', err.message)
end

if isa(alpha, 'skpParameter') && alpha.state == paramState.onInnerBdry
    fprintf('On inner boundary %d. Hejhal not needed.\n', alpha.ison)
    return
end

ja = find(abs(alpha - dv) <= qv);
fprintf('Parameter is inside circle %d.\n', ja);


%%
% Step 2: Define a "pair parameter" beta in the domain, and find logXHat
% from this parameter.

beta = theta(ja, 1/conj(alpha)); % alpha = theta(ja, 1/conj(beta))!
wb = skprime(beta, dv, qv);

logXhatBeta = @(z) wb.evalLogXhat(z);


%%
% Step 3: Use this to define logXhat for alpha via Hejhal's
% identity.

vj = wb.vjFuns{ja};
% logXhatAlpha = @(z) -4i*pi*(vj(1/conj(alpha)) - vj(zeta) + vj.taujj/2) ...
%     + 2*log(qv(ja)*conj(alpha)/conj(alpha - dv(ja))...
%             *(z - alpha)./(z - theta(ja, 1/conj(alpha)))) ...
%     + logXhatBeta(z);



































