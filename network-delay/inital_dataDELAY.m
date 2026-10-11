function[U0,Z0]=inital_dataDELAY(N,Nx,Nrho,L,h)
size=Nx-1;
dx2=L/(Nx-2);
drho=1/(Nrho-1);
x=0:dx2:L;
rho=0:drho:1;
    for i=1:N;
        U0((i-1)*size+1:i*size)=(1-cos(2*pi*x/L));
        Z0(1:size,(i-1)*Nrho+1:i*Nrho)=(1-cos(2*pi*x/L))'*cos(2*pi*rho*h(i));
    end
    U0=U0';
    Z0=Z0;
end

