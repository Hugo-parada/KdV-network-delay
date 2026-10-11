clear all
clc
% this code solve a kdv on network with internal delay

N=4;
T=20;
L=2; %we choose the same length for all branches
alpha=3;
Nx=100; %Nx 100 works well take Nx such that Nx/10,5,4,2 are integers
Nt=100;  %Nt=50
dx=L/Nx; %space discretization
dt=T/Nt; %time discretization

%%For this implementation the set omega will be seting as (0,L/2),(0,L/4)
% (0,L/2) and (0,L/4)

h=[2 1 2 1];  %delay for each branch 


%caso 1
a=[ 0 0 0 0];  %feedback for each branch
b=[ 0 0 0 0];  %delayed feedback for each branch 
xi=[0 0 0 0];


%caso 2
%a=[ 2 1 1 1];  %feedback for each branch
%b=[ 0 0 0 0];  %delayed feedback for each branch 
%xi=[0 0 0 0];


%caso 3
%a=[ 1 1 1 1];  %feedback for each branch
%b=[ 1 1 1 1];  %delayed feedback for each branch 
%xi=[1 1 1 1];


%caso 4
%a=[ 5 4 3 0];  %feedback for each branch
%b=[ 2 1 2 0.5];  %delayed feedback for each branch 
%xi=[3 2 3 1];


%caso 5
%a=[ 0 5 3 0];  %feedback for each branch
%b=[8 2 1 2];  %delayed feedback for each branch 
%xi=[8 3 2 3];



[Aj,Bj]=omega(N,Nx,a,b);

%%

%rho is moving bewteen 0,1
Nrho=5;
drho=1/Nrho;





[AN,Dplus,Dz,vec]=MatDELAY(N,Nx,Nrho,L,dt,alpha);
[U0,Z0]=inital_dataDELAY(N,Nx,Nrho,L,h);
fullU0=fullsol(U0,N,Nx);
U=fullU0;
%U=U0;
vecz=zeros(Nrho,1);
vecz(1)=1;
Z=Z0;


Aint=int16.empty(Nx-1,0);

Eu=0;
Ez=0;
EU=Eu;
EZ=Ez;
for n=1:Nt;
    Zn= double.empty(size(Aint));
    %Zrho1=double.empty(size(int16.empty(Nx-1,0)));
    Zrho1=double.empty(size(int16.empty(0,1)));
    
    EU=[EU dx*norm(U0)^2];
    
    for j=1:N;
        
        Ze=Z0(1:Nx-1,(j-1)*Nrho+1:j*Nrho);
        Ze=Ze(:);
        Ez=Ez+h(j)*xi(j)*dx*drho*norm(Ze)^2;
    
        Zbc=U0((j-1)*(Nx-1)+1:j*(Nx-1));
        Zi=(h(j)*eye(Nrho)+dt*Dz)\(   (h(j)/dt)*Z0(1:Nx-1,(j-1)*Nrho+1:j*Nrho)'+(dt/drho)*vecz*Zbc'   );
        aux=Zi';
        Zn=[Zn, aux];
        Zrho1=[Zrho1;aux(1:end,Nrho)];
    end
    Z0=Zn;
    Z=[Z;Z0];
    EZ=[EZ Ez];
    Ez=0;
    
    
    NL=nonlinearpart(U0,Dplus,N);
    
    %Without delay without feedback
    %Ui=AN\((1/dt)*(U0.*vec)+NL);  
    
    %Without delay with feedback
    %Ui=(AN+diag(Aj))\((1/dt)*(U0.*vec)+NL);  
    
    %delay without feedback
    %Ui=AN\((1/dt)*(U0.*vec)+NL-diag(Bj)*Zrho1);  
    
    %delay  feedback
    %Ui=(AN+diag(Aj))\((1/dt)*(U0.*vec)+NL-diag(Bj)*Zrho1);  %ESA 
    
    %% esta es la prueba con el dt por fuera
    %Ui=(AN+dt*diag(Aj))\(  (U0.*vec)+dt*NL-dt*diag(Bj)*Zrho1 );
    
    % otra prueba
    %Ui=(AN)\(  (U0.*vec)+dt*NL-dt*diag(Bj)*Zrho1-dt*diag(Aj)*U0 );
    %Ui=AN\((1/dt)*(U0.*vec)+NL-diag(Bj)*Zrho1-diag(Aj)*U0); 
    
    % prueba full explicit
    %Ui=(AN*U0-dt*diag(Aj)*U0-dt*diag(Bj)*Zrho1-diag(Aj)*U0-)
    
    U0=Ui;
    fullU0=fullsol(U0,N,Nx);
    U=[U fullU0];
    %U=[U U0];
end 
EU=EU(2:end);
EZ=EZ(2:end);
E=EU+EZ;

[x,t]=meshgrid(0:dx:L,0:dt:T);
surf(U)
%surf(t,x,U')E=norm(U(:,1));




t2=0:T/(Nt-1):T;
%loglog(t2,E)

%plot(t2,log(E))