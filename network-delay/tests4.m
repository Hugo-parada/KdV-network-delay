%% Figure 5: numerical decay and the Theorem 1.1 bound
clear; clc;
T=2; Nt=100; L=2; alpha=3; N=4;
h=[0.1,0.2,0.3,0.4];
a=[2,4,6,8]; b=[0.5,1.5,2.5,3.5]; xi=[1,2,3,4];

theta=0.1;
[E,t]=plotkdvdnet(a,b,xi,h,T,Nt,theta);
epsilon=1.01*sqrt(E(1));
[gamma,C,info]=theorem11_bound(a,b,xi,h,L,alpha,N,epsilon);

theoretical=log(C*E(1))-2*gamma*t;
valid=isfinite(E) & E>0;
p=polyfit(t(valid),log(E(valid)),1);
numerical_fit=polyval(p,t);

figure;
plot(t,log(E),'LineWidth',1.2); hold on;
plot(t,theoretical,'--','LineWidth',1.2);
plot(t,numerical_fit,':','LineWidth',1.2); hold off;
legend('numerical energy','Theorem 1.1 bound','linear fit','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;

fprintf('epsilon = %.8g (epsilon_max = %.8g)\n',epsilon,info.epsilon_max);
fprintf('mu1 = %.8g (upper bound %.8g)\n',info.mu1,info.mu1_max);
fprintf('mu2 = %.8g (upper bound %.8g)\n',info.mu2,info.mu2_max);
fprintf('gamma = %.8g\n',gamma);
fprintf('C = %.8g\n',C);
fprintf('linear fit = %.8g + (%.8g)t\n',p(2),p(1));
