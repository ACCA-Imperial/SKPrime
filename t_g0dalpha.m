%% Testing G0 alpha derivative.
clear

dv = [
    -0.2517+0.3129i
    0.2307-0.4667i
    ];
qv = [
    0.2377
    0.1557
    ];
D = skpDomain(dv, qv);

alpha = 1e6*exp(1i*pi/4);


%%

dag0 = greensC0da(alpha, D);

plotfd(dag0, D, 'unit')
