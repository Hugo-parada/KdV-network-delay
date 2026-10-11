% Nonlinear KdV on a star-shaped network without delay.
clear; clc;

N = 3;
L = 5;
alpha = 3;
T = 20;
Nt = 400;
Nx = 100;
theta = 1;

[t,E,out] = solve_kdv_network(N,L,alpha,T,Nt,Nx,theta);

figure;
for j = 1:N
    subplot(N,1,j);
    surf(out.x,t,squeeze(out.U(:,j,:))');
    xlabel('x'); ylabel('t'); zlabel(sprintf('u_%d',j));
    shading interp;
end

figure;
plot(t,log(E),'LineWidth',1.2);
xlabel('time'); ylabel('ln(E(t))'); grid on;
