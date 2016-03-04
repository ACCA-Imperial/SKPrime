%% G0 alpha FD derivative.
clear

D = dat_test_data(3);
alpha = 0.5;

h = 1e-6;

g0 = greensC0(alpha, D);


%%

da = alpha + h*[1, -1, 1i, -1i]/2;

j = 1;
drt = greensCj(da(1), j, D);
dlf = greensCj(da(2), j, drt);
dup = greensCj(da(3), j, drt);
ddn = greensCj(da(4), j, drt);

dgj = @(z) (drt.hat(z) - dlf.hat(z) + dup.hat(z) - ddn.hat(z))/h;


%%

% plotfd(dg0, D, 'unit')


%%

zb = boundaryPts(D, 50);
zb = zb(:,1:3);

imag(dgj(zb)) + real(1./(zb - alpha))/(2*pi)
