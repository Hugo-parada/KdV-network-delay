%% First simulations


clear all
clc


T=5;
Nt=100;  

h=[1 0.5 0.5 1];  %delay for each branch


%caso 1

a_1=[ 0 0 0 0];  %feedback for each branch
b_1=[ 0 0 0 0];  %delayed feedback for each branch 
xi_1=[0 0 0 0];


%caso 2

a_2=[ 10 10 10 10];  %feedback for each branch
b_2=[ 0 0 0 0];  %delayed feedback for each branch 
xi_2=[0 0 0 0];


%caso 3

a_3=[ 1 1 1 1];  %feedback for each branch
b_3=[ 1 1 1 1];  %delayed feedback for each branch 
xi_3=[1 1 1 1];


%caso 4

a_4=[ 1 1 1 0];  %feedback for each branch
b_4=[ 1 1 1 0.5];  %delayed feedback for each branch 
xi_4=[3 2 3 1];


%caso 5


a_5=[ 0 0 0 0];  %feedback for each branch
b_5=[10 10 10 10];  %delayed feedback for each branch 
xi_5=b_5;

%xi_5=[15 15 15 15];
big=1;

[E_1,t_1]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,big);
plot(t_1,log(E_1))

[E_2,t_2]=plotkdvdnet(a_2,b_2,xi_2,h,T,Nt,big);
%plot(t_2,log(E_2))

[E_3,t_3]=plotkdvdnet(a_3,b_3,xi_3,h,T,Nt,big);
%plot(t_3,log(E_3))

[E_4,t_4]=plotkdvdnet(a_4,b_4,xi_4,h,T,Nt,big);
%plot(t_4,log(E_4))

[E_5,t_5]=plotkdvdnet(a_5,b_5,xi_5,h,T,Nt,big);
%plot(t_5,log(E_5))






plot(t_1,log(E_1),'-',t_1,log(E_2),'--',t_1,log(E_3),'-.',t_1,log(E_4),':',t_1,log(E_5),'.')
legend('a_{j}=b_{j}=0','a_{j}=10, b_{j}=0','a_{j}=b_{j}=1','a_{1}=a_{2}=a_{3}=1, a_{4}=0, b_{1}=b_{2}=b_{3}=1, b_{4}=0.5,','a_{j}=0, b_{j}=10')
xlabel('time ') 
ylabel('ln(E)')
