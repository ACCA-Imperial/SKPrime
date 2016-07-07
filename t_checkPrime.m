%% Prime test script.

clear

L = 6;


%%

td = skpUnitTest.domainSimple3;
dv = td.dv;
qv = td.qv;

wp = skprod(dv, qv, L);


%%

alpha = td.parameterNearCirc;

tp = [td.testPointInside;
    td.testPointOutside];

w = skprime(alpha, dv, qv);

error = wp(tp, alpha) - w(tp);

disp(error)
