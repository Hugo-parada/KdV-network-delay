function[U0,Z0]=inital_dataDELAY(N,Nx,Nrho,L,h)

size=Nx-1;
dx=L/Nx;
drho=1/Nrho;

% the stored KdV points are x_0,...,x_{Nx-2}
x=(0:size-1)*dx;
% rho=0 is the inflow z(t,0,x)=u(t,x), so Z stores rho_1,...,rho_Nrho
rho=(1:Nrho)*drho;

U0=zeros(N*size,1);
Z0=zeros(size,N*Nrho);

for j=1:N;
    U0((j-1)*size+1:j*size)=(1-cos(2*pi*x/L))';
    Z0(:,(j-1)*Nrho+1:j*Nrho)=...
        (1-cos(2*pi*x/L))'*cos(2*pi*rho*h(j));
end
end
