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

alpha0 = (0.53528 + 0.17493i);
% alpha = dv(1) + qv(1);

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);

alpha = theta(j, 1/conj(alpha0)); % <-------- the parameter of concern


%%
% Step 1: Identify which circle the parameter is in. Is it a "first level
% reflection" from the fundamental domain?

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

% thj = @(z) dv(ja) + qv(ja)^2*z./(1 - conj(dv(ja))*z);
beta = theta(ja, 1/conj(alpha));
if ~(isin(D, beta) || isin(D, 1/conj(beta)))
    fprintf('Parameter is more than one Schottky group level deep.\n')
    return
end


%%
% Step 2: Get the Green's function for the hole and use it to get the
% desired prime function.

wb = skprime(beta, dv, qv);
Gj = greensCj(beta, ja, wb);

w = @(z) wb(z)./exp(2i*pi*Gj(z))/abs(beta - dv(ja));



































