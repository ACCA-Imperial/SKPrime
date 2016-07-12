%% How to test the Gj?
clear

j = 1;
L = 6;


%%

td = skpUnitTest.domainSimple3;
D = skpDomain(td);
% [dv, qv] = domainData(D);
dj = D.dv(j);
qj = D.qv(j);

alpha = td.parameter('origin');
tp = td.testPoint('inside');

thj = @(z) D.theta(j, z);


%%

gj = greensCj(alpha, j, D);
egj = @(z) exp(2i*pi*gj(z));
ehgj = @(z) exp(2i*pi*gj.hat(z));

wtop = skprime(alpha, gj);
wbot = skprime(thj(1/conj(alpha)), wtop);

refj = @(z) wtop(z)./wbot(z)*qj/abs(alpha - dj);
refh = @(z) refj(z)./(z - alpha).*(z - thj(1/conj(alpha)));

% err = refj(tp) - egj(tp);
err = refh(tp) - ehgj(tp);
disp(err)
