% Example: nonlinear KdV network with localized internal delay.
clear; clc;

T = 5;
Nt = 100;
h = [1,0.5,1,1];
a = [1,1,1,1];
b = [1,1,1,1];
xi = [1,1,1,1];
theta = 1;

[E,t,out] = plotkdvdnet(a,b,xi,h,T,Nt,theta);

figure;
plot(t,log(E),'LineWidth',1.2);
xlabel('time'); ylabel('ln(E(t))'); grid on;

figure;
for j = 1:4
    subplot(4,1,j);
    surf(out.x,t,squeeze(out.U(:,j,:))');
    xlabel('x'); ylabel('t'); zlabel(sprintf('u_%d',j));
    shading interp;
end
