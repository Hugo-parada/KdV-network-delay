%% Figure 2: feedback configurations
clear; clc;
T = 5; Nt = 100;
h = [1,0.5,1,1];

a_1=[0,0,0,0]; b_1=[0,0,0,0]; xi_1=[0,0,0,0];
a_2=[10,10,10,10]; b_2=[0,0,0,0]; xi_2=[0,0,0,0];
a_3=[1,1,1,1]; b_3=[1,1,1,1]; xi_3=[1,1,1,1];
a_4=[1,1,1,0]; b_4=[1,1,1,0.5]; xi_4=[3,2,3,1];
a_5=[0,0,0,0]; b_5=[10,10,10,10]; xi_5=b_5;
theta=1;

[E_1,t]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,theta);
[E_2,~]=plotkdvdnet(a_2,b_2,xi_2,h,T,Nt,theta);
[E_3,~]=plotkdvdnet(a_3,b_3,xi_3,h,T,Nt,theta);
[E_4,~]=plotkdvdnet(a_4,b_4,xi_4,h,T,Nt,theta);
[E_5,~]=plotkdvdnet(a_5,b_5,xi_5,h,T,Nt,theta);

figure;
plot(t,log(E_1),'-',t,log(E_2),'--',t,log(E_3),'-.', ...
     t,log(E_4),':',t,log(E_5),'.-','LineWidth',1.1);
legend('a_j=b_j=0','a_j=10, b_j=0','a_j=b_j=1', ...
       'a_{1,2,3}=1, a_4=0; b_{1,2,3}=1, b_4=0.5', ...
       'a_j=0, b_j=10','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
