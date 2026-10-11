%% Figure 3: influence of delayed feedback and delay sizes
clear; clc;
T=10; Nt=200; theta=1;

a_1=[0,0,0,0]; b_1=[20,20,20,20]; xi_1=b_1; h_1=[0,0,0.5,1];
a_2=[0,0,0,0]; b_2=[20,20,20,20]; xi_2=b_2; h_2=[1,1,1,1];
a_3=[0,0,0,0]; b_3=[25,25,25,25]; xi_3=b_3; h_3=[1,1,1,1];

[E_1,t]=plotkdvdnet(a_1,b_1,xi_1,h_1,T,Nt,theta);
[E_2,~]=plotkdvdnet(a_2,b_2,xi_2,h_2,T,Nt,theta);
[E_3,~]=plotkdvdnet(a_3,b_3,xi_3,h_3,T,Nt,theta);

figure;
plot(t,log(E_1),'-',t,log(E_2),'--',t,log(E_3),'-.','LineWidth',1.1);
legend('h=(0,0,0.5,1), b_j=20','h_j=1, b_j=20', ...
       'h_j=1, b_j=25','Location','best');
xlabel('time'); ylabel('ln(E(t))'); grid on;
