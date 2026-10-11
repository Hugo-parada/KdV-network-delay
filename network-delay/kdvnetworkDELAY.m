% Example of the delayed KdV network
clear;clc;

a=[2 4 6 8];
b=[0.5 1.5 2.5 3.5];
xi=[1 2 3 4];
h=[0.1 0.2 0.3 0.4];
T=10;
Nt=500;
big=0.5;

[E,t]=plotkdvdnet(a,b,xi,h,T,Nt,big);
plot(t,log(E),'LineWidth',1.2)
xlabel('time')
ylabel('ln(E(t))')
grid on
