%% G0 product checks.
clear

L = 6;


%%

dv = [0; 0.52621+0.23899i];
qv = [0.25; 0.14095];
D = skpDomain(dv, qv);

parameterInside = -0.50729+0.29388i;
testPointInside = 0.27638-0.55977i;

alpha = parameterInside;

zp = testPointInside;
zb = boundaryPts(D, 5);


%%

wp = skprod(dv, qv, L);

logprat = @(z,a) log(wp(z,a)./wp(z, 1/conj(a)))/2i/pi;
if alpha ~= 0
    g0p = @(z,a) logprat(z, a) - log(abs(a))/2i/pi;
    g0hp = @(z,a) g0p(z, a) - log((z - a)./(z - 1/conj(a)))/2i/pi;
else
    g0p = @(z,a) logprat(z, a);
    g0hp = @(z,a) g0p(z, a) - log(z)/2i/pi;
end


%%

g0 = greensC0(alpha, D);

err = exp(2i*pi*g0p(zb, alpha)) - exp(2i*pi*g0(zb));
disp(abs(err))

errh = exp(2i*pi*g0hp(zb, alpha)) - exp(2i*pi*g0.hat(zb));
disp(abs(errh))

































