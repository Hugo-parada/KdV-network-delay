function[AN,Dplus,Dz,vec]=MatDELAY(N,Nx,Nrho,L,dt,alpha)

size=Nx-1;
dx=L/Nx;
drho=1/Nrho;

beta1=(-1/dx^3-1/(2*dx));
beta2=3/dx^3;
beta3=(-3/dx^3+1/(2*dx));
beta4=1/dx^3;

A=diag(beta2*ones(1,size))+diag(beta3*ones(1,size-1),1)...
  +diag(beta1*ones(1,size-1),-1)+diag(beta4*ones(1,size-2),2);
A(1,1:3)=0;

AN=zeros(N*size);
vec=ones(N*size,1);
Dplus=zeros(N*size);

for j=1:N;
    pos=(j-1)*size+1:j*size;
    AN(pos,pos)=A;

    % Forward difference within each branch.
    Dpj=(1/dx)*(diag(-ones(1,size))+diag(ones(1,size-1),1));
    Dpj(1,:)=0;
    Dplus(pos,pos)=Dpj;

    vec((j-1)*size+1)=0;
end

% Kirchhoff condition at the common node:
% sum_j u_j,xx(0)+alpha*u(0)=-(N/3)u(0)^2.
% Keeping the N copies of the node gives N/dx^2+alpha after continuity.
AN(1,:)=0;
for j=1:N;
    p=(j-1)*size;
    AN(1,p+1)=AN(1,p+1)+1/dx^2;
    AN(1,p+2)=AN(1,p+2)-2/dx^2;
    AN(1,p+3)=AN(1,p+3)+1/dx^2;
end
AN(1,1)=AN(1,1)+alpha;

% continuity at the node
for j=2:N;
    row=(j-1)*size+1;
    AN(row,:)=0;
    AN(row,row)=1;
    AN(row,1)=-1;
end

AN=dt*AN+eye(N*size)*diag(vec);
Dz=(1/drho)*(diag(ones(1,Nrho))+diag(-ones(1,Nrho-1),-1));
end
