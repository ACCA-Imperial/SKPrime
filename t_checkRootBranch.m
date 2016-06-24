%% Check root branch.
clear

% j = 2;
L = 6;
tol = 5e-3;


%%

D = fd_test_data(5);
[dv, qv, m] = domainData(D);

zp = [
    -0.46531-0.20292i
    -0.12945+0.43382i
    0.58426-0.027988i];
zp = [zp; 1./conj(zp)];

% dv = [
%   -0.2517+0.3129i
%    0.2307-0.4667i];
% qv = [
%   0.2377
%   0.1557];
% D = skpDomain(dv, qv);
% 
% m = numel(dv);
%
% zp = [0.5; 2];

theta = @(j,z) dv(j) + qv(j)^2*z./(1 - conj(dv(j))*z);



%%

fprintf('Setting up skprod ...\n')
wref = skprod(dv, qv, L);
fprintf('\b done\n')


%%
% Place "generator" parameter at points around each circle looking for
% branch cut problems.
%
% TODO: Use test points inside and outside unit region.

np = 40;
for j = 1:m
    dj = dv(j);
    qj = qv(j);
    
    for ij = 1:m
        av = 2*pi*(0:np-1)'/(np - 1);
        av = dv(ij) + (qv(ij) + 0.1)*exp(1i*av);
        
        parfor i = 1:np
            fprintf('%d,%d,%d\n', j, ij, i)

            alpha = av(i);
            beta = 1/conj(alpha);
            thja = theta(j, beta);
            
            wb = skprime(beta, D);
            vj = wb.vjFuns{j};
            taujj = vj.taujj;
            rootHj = @(z) -exp(-2i*pi*(vj(beta) - vj(z) + taujj/2)) ...
                *qj/(1 - conj(dj)*beta);
            wja = @(z) wb(z).*rootHj(z);
            
            aerr = abs(wref(zp, thja) - wja(zp)); %#ok<PFBNS>
            if any(aerr > tol)
                je = find(aerr > tol, 1);
                fprintf(['(j=%d) point idx %d near circ %d for angle ' ...
                    '%.4f*pi, abserr=%.4g.\n'], ...
                    j, je, ij, (i - 1)/(np - 1)*2, aerr(je))
            end
        end
    end
end




























