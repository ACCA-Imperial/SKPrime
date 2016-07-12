%% How to test the Gj?
clear

j = 1;
L = 6;


%%

td = skpUnitTest.domainSimple3;
D = skpDomain(td);
[dv, qv] = domainData(D);

alpha = td.parameter('inside');
tp = td.testPoint('inside');

thj = @(z) D.theta(j, z);


%%

gj = greensCj(alpha, j, D);
egj = @(z) exp(2i*pi*gj(z));
ehgj = @(z) exp(2i*pi*gj.hat(z));

wtop = skprime(alpha, gj);
wbot = skprime(thj(1/conj(alpha)), wtop);

ref1 = @(z) wtop(z)./wbot(z)*qv(j)/abs(alpha - dv(j));
refh = @(z) ref1(z)...
     ./(z - alpha).*(z - thj(1/conj(alpha)))/qv(j)*abs(alpha - dv(j));

err = refh(tp) - ehgj(tp);
disp(err)
