function[E,t2]=plotkdvdnet(a,b,xi,h,T,Nt,big)
% this code solves the nonlinear KdV equation on a star network with
% localized instantaneous feedback and internal delay

N=4;
L=2;
alpha=3;
Nx=100;
Nrho=100;

dx=L/Nx;
dt=T/Nt;
drho=1/Nrho;

[Aj,Bj,O]=omega(N,Nx,a,b);
[AN,Dplus,Dz,vec]=MatDELAY(N,Nx,Nrho,L,dt,alpha);

vec2=vec;
vec2(1)=1;
vecz=zeros(Nrho,1);
vecz(1)=1;

[U0,Z0]=initial_dataDELAY(N,Nx,Nrho,L,h);
U0=big*U0;
Z0=big*Z0;

size=Nx-1;
t2=0:dt:T;
E=zeros(1,Nt+1);

% initial energy
Eu=dx*norm(U0)^2;
Ez=0;
for j=1:N;
    Ze=Z0(:,(j-1)*Nrho+1:j*Nrho);
    Oj=O((j-1)*size+1:j*size)';
    Ez=Ez+h(j)*xi(j)*dx*drho*sum(sum((Oj*ones(1,Nrho)).*(Ze.^2)));
end
E(1)=Eu+Ez;

for n=1:Nt;
    Zn=zeros(size,N*Nrho);
    Zrho1=zeros(N*size,1);

    % transport equation h_j z_t+z_rho=0 on each branch
    for j=1:N;
        Zbc=U0((j-1)*size+1:j*size);
        Zold=Z0(:,(j-1)*Nrho+1:j*Nrho);
        Zi=(h(j)*eye(Nrho)+dt*Dz)\...
            (h(j)*Zold'+(dt/drho)*vecz*Zbc');
        Zi=Zi';
        Zn(:,(j-1)*Nrho+1:j*Nrho)=Zi;
        Zrho1((j-1)*size+1:j*size)=Zi(:,Nrho);
    end
    Z0=Zn;

    NL=nonlinearpart(U0,Dplus,N);

    % implicit linear KdV part, explicit nonlinear part
    Ui=(AN+dt*diag(Aj)*diag(vec))\...
       ((U0-dt*diag(Bj)*Zrho1).*vec-dt*NL.*vec2);
    U0=Ui;

    Eu=dx*norm(U0)^2;
    Ez=0;
    for j=1:N;
        Ze=Z0(:,(j-1)*Nrho+1:j*Nrho);
        Oj=O((j-1)*size+1:j*size)';
        Ez=Ez+h(j)*xi(j)*dx*drho*sum(sum((Oj*ones(1,Nrho)).*(Ze.^2)));
    end
    E(n+1)=Eu+Ez;
end
end
