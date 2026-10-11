%% Second simulations


clear all
clc



T=10;
Nt=200; 
big=1;

%caso 1
a_1=[ 0 0 0 0];  %feedback for each branch
b_1=[20 20 20 20];  %delayed feedback for each branch 
xi_1=b_1;
h_1=[0 0 0.5 1];  %delay for each branch


%caso 2
a_2=[ 0 0 0 0];  %feedback for each branch
b_2=[20 20 20 20];  %delayed feedback for each branch 
xi_2=b_2;
h_2=[1 1 1 1];  %delay for each branch


%caso 3
a_3=[ 0 0 0 0];  %feedback for each branch
b_3=[25 25 25 25];  %delayed feedback for each branch %25 explote
xi_3=b_1;
h_3=[1 1 1 1];  %delay for each branch



[E_1,t_1]=plotkdvdnet(a_1,b_1,xi_1,h_1,T,Nt,big);
plot(t_1,log(E_1))

[E_2,t_2]=plotkdvdnet(a_2,b_2,xi_2,h_2,T,Nt,big);
%plot(t_2,log(E_2))

[E_3,t_3]=plotkdvdnet(a_3,b_3,xi_3,h_3,T,Nt,big);
%plot(t_3,log(E_3))

zero=0*t_1;




plot(t_1,log(E_1),'-',t_1,log(E_2),'--',t_1,log(E_3),'-.')
legend('h_{1}=h_{2}=0, h_{3}=0.5, h_{4}=1, b_{j}=20','h_{j}=1, b_{j}=20','h_{j}=1, b_{j}=25')
xlabel('time ') 
ylabel('ln(E)')
