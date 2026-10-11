clear all
clc
% this code solve a kdv on network
N=3;
T=50;
L=5; %we choose the same length for all branches
alpha=3;
Nx=10; 
Nt=5;
dx=L/Nx; %space discretization
dt=T/Nt; %time discretization

[AN,Dplus,vec]=Mat(N,Nx,L,dt,alpha);
U0=inital_data(N,Nx,L);
fullU0=fullsol(U0,N,Nx);
U=fullU0;
%U=U0;
for n=1:Nt;
    NL=nonlinearpart(U0,Dplus,N);
    Ui=AN\((1/dt)*(U0.*vec)+NL);
    U0=Ui;
    fullU0=fullsol(U0,N,Nx);
    U=[U fullU0];
    %U=[U U0];
end

%[x,t]=meshgrid(0:dx:L,0:dt:T);
surf(U)
%surf(t,x,U')
  
        