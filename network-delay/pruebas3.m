%% Figure 4: normalized energy for different initial amplitudes
clear; clc;
T=2; Nt=100;
h=[0.1,0.2,0.3,0.4];
a=[2,4,6,8]; b=[0.5,1.5,2.5,3.5]; xi=[1,2,3,4];
theta_values=[0.5,2,5];

[E_1,t]=plotkdvdnet(a,b,xi,h,T,Nt,theta_values(1));
[E_2,~]=plotkdvdnet(a,b,xi,h,T,Nt,theta_values(2));
[E_3,~]=plotkdvdnet(a,b,xi,h,T,Nt,theta_values(3));

figure;
plot(t,log(E_1/E_1(1)),'-',t,log(E_2/E_2(1)),'--', ...
     t,log(E_3/E_3(1)),'-.','LineWidth',1.1);
legend('\theta=0.5','\theta=2','\theta=5','Location','best');
xlabel('time'); ylabel('ln(E(t)/E(0))'); grid on;
