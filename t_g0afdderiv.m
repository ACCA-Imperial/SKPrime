%% Gj alpha FD derivative.
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
alpha = 0.5;

h = 1e-6;


%%

j = 0;
if j == 0
    gfun = @greensC0;
else
    gfun = @(a, D) greensCj(a, j, D);
end

gj = gfun(alpha, D);

da = alpha + h*[-1, 1]/2;
gh1 = gfun(da(1), gj);
gh2 = gfun(da(2), gj);

% Alpha derivative of Gj.
dagj = @(z) (gh2(z) - gh1(z))/h;

% Alpha derivative of Gjhat.
dagjh = @(z) (gh2.hat(z) - gh1.hat(z))/2/h;


%%

[zb, t] = boundaryPts(D, 200);
zb = zb(:,1:3);

clf
% plot(t, abs(real(1./(zb - alpha))/2/pi - ...
%     imag(dagj(zb)) - imag(dagjh(zb))) )
plot(t, imag(dagjh(zb)), '-')
hold on
% set(gca, 'colororderindex', 1)
plot(t, real(1./(zb - alpha))/2/pi, '--')
hold off
xlim([0, 2*pi])


%%

return

% Points on D_zeta circles.
zb = boundaryPts(D, 25);
zb = zb(:,1:3);

% Constant.
disp(imag(gj(zb)))

% Also constant.
% disp(imag(dagj(zb)))

% NOT CONSTANT !?!?
disp(imag(dagjh(zb)) + real(1./(zb - alpha))/2/pi)
