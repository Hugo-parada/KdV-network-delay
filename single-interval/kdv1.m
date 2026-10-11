% Nonlinear KdV on one bounded interval.
clear; clc;

L = 5;
alpha = 3;
T = 10;
Nt = 200;
Nx = 100;
theta = 1;

[t,E,out] = solve_single_interval(L,alpha,T,Nt,Nx,theta);

figure;
surf(out.x,t,squeeze(out.U(:,1,:))');
xlabel('x'); ylabel('t'); zlabel('u(t,x)');
shading interp;

figure;
plot(t,log(E),'LineWidth',1.2);
xlabel('time'); ylabel('ln(E(t))'); grid on;
