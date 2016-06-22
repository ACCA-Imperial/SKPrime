%% Check points under theta (again).
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

av = [
    -0.34636+0.2449i
    0.38134+0.0069971i
    ];
av = [av; 1./conj(av)];

thj = @(z) D.theta(j, z);


%%

mask = isin(D, av) | isin(D, 1./conj(av));
for i = find(mask(:))'
    fprintf(' av(%d) is in FD.\n', i)
end

mask = isin(D, thj(1./conj(av)));
for i = find(mask(:))'
    fprintf(' av(%d) is 1st level in disk C_j\n', i)
end

mask = isin(D, thj(av));
for i = find(mask(:))'
    fprintf(' av(%d) is 1st level in disk C_j''\n', i)
end
