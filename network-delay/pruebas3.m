clear all
clc


T=2;
Nt=100;  

h=[0.1 0.2 0.3 0.4];  %delay for each branch


%caso 1

a_1=[ 2 4 6 8];  %feedback for each branch
b_1=[ 0.5 1.5 2.5 3.5];  %delayed feedback for each branch 
xi_1=[1 2 3 4];

big_1=0.5;
big_2=2;
big_3=5;




[E_1,t_1]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,big_1);

[E_2,t_2]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,big_2);

[E_3,t_3]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,big_3);







plot(t_1,log(E_1/E_1(1)),'-',t_1,log(E_2/E_2(1)),'--',t_1,log(E_3/E_3(1)),'-.')
legend('\theta=0.5','\theta=2','\theta=5')
xlabel('Time ') 
ylabel('ln(E)')
