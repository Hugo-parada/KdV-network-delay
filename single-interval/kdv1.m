clear all
clc
% this code solve a kdv 
T=10;
L=5; %we choose the same length for all branches
alpha=3;
Nx=100; 
Nt=50;
dx=L/Nx; %space discretization
dt=T/Nt; %time discretization

%%
[A,Dplus]=Mat1(Nx,L,dt,alpha);
U0=inital_data1(Nx,L);
fullU0=[U0;0;0];
U=fullU0;
vec=ones((Nx-1),1);
vec(1)=0;
for n=1:Nt;
    NL=nonlinearpart1(U0,Dplus);
    Ui=A\((1/dt)*(U0.*vec)+NL);
    U0=Ui;
    fullU0=[U0;0;0];
    U=[U fullU0];
end

[x,t]=meshgrid(0:dx:L,0:dt:T);
%surf(U)
surf(t,x,U')

