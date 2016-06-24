%% Check points under theta (again).
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

av = [
    -0.34636+0.2449i
    0.38134+0.0069971i
    ];
av = [av; 1./conj(av)];

thj = @(z) D.theta(j, z);


%%

for i = 1:numel(av)
    ap = skpParameter(av(i), D);
    fprintf(' av(%d) is in %s', i, char(ap.state))
    if ap.state == paramState.innerDisk ...
            || ap.state == paramState.outerDisk
        fprintf(' (#%d)', ap.indisk)
    end
    fprintf('.\n')
end
