%% Check root branch.
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

m = numel(dv);

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);

zp = 0.5;
wref = skprod(dv, qv, 6);


%%
% Place "generator" parameter at points around each circle looking for
% branch cut problems.
%
% TODO: Use test points inside and outside unit region.

np = 35;
for ij = 1:m
    av = 2*pi*(0:np-1)'/(np - 1);
    av = dv(ij) + (qv(ij) + 0.1)*exp(1i*av);
    
    for i = 1:np
        alpha = av(i);
        beta = 1/conj(alpha);        
        thja = theta(j, beta);
        
        wb = skprime(beta, dv, qv);
        vj = wb.vjFuns{j};
        taujj = vj.taujj;
        rootHj = @(z) -exp(-2i*pi*(vj(beta) - vj(z) + taujj/2)) ...
            *qv(j)/(1 - conj(dv(j))*beta);
        wja = @(z) wb(z).*rootHj(z);
        
        if abs(wref(zp, thja) - wja(zp)) > 1e-4
            fprintf('\nAlpha near circle %d at angle %.4f/pi looks to be a problem.\n', ...
                ij, (i - 1)/(np - 1)*2)
        else
            fprintf('.')
        end
    end
end
fprintf('\n')






























