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

% alpha = 0.53528 + 0.17493i;
alpha = dv(1) + qv(1);

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);

aj = theta(j, 1/conj(alpha)); % <-------- the parameter of concern


%%
% Step 1: Identify which circle the parameter is in.

try
    aj = skpParameter(aj, D);
catch err
    fprintf('%s\n', err.message)
end

if isa(aj, 'skpParameter') && aj.state == paramState.onInnerBdry
    fprintf('On inner boundary %d. Hejhal not needed.\n', aj.ison)
    return
end

ja = find(abs(aj - dv) <= qv);
fprintf('Parameter is inside circle %d.\n', ja);
