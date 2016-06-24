%% Try to extend v_j function into holes.
clear


%%

dv = [
  -0.2517+0.3129i
   0.2307-0.4667i];
qv = [
  0.2377
  0.1557];

D = skpDomain(dv, qv);
[~, ~, m, di, qi] = domainData(D);

j = 1;

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);


%%

vj = vjFirstKind(j, D);
tjj = vj.taujj;


%%

res = 800;
scale = 1.5;
zg = repmat(linspace(-scale, scale, res), res, 1);
zg = complex(zg, zg');

val = vj(zg);
% val = complex(nan(size(zg)));


%%
% Points in D_zeta.

% mask = isin(D, zg);
% val(mask) = vj(zg(mask));
% 
% 
% %%
% % Points theta_j(1/conj(z)) in D_zeta. (1st reflection into C_j.)
% 
% done = mask;
% mask(mask) = false;
% mask(~done) = isin(D, theta(j, 1./conj(zg(~done))));
% val(mask) = conj(vj(theta(j, 1./conj(zg(mask)))) - tjj);
% 
% 
% %%
% % Points 1/z in D_zeta.
% 
% done = done | mask;
% mask(mask) = false;
% mask(~done) = isin(D, 1./conj(zg(~done)));
% val(mask) = conj(vj(1./conj(zg(mask))));
% 
% 
% %%
% % Points theta_j(z) in D_zeta. (1st reflection into C_j'.)
% 
% done = done | mask;
% mask(mask) = false;
% mask(~done) = isin(D, theta(j, zg(~done)));
% val(mask) = vj(theta(j, zg(mask))) - tjj;


%%

figure(1), clf
ezPhase(zg, val, 'c')
hold on
plot(D)
for i = 1:m
    plot(circle(di(i), qi(i)))
end
hold off
axis(scale*[-1, 1, -1, 1])


































