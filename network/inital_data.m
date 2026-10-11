function[U0]=inital_data(N,Nx,L)
size=Nx-1;
dx2=L/(Nx-2);
x=0:dx2:L;
    for i=1:N;
        U0((i-1)*size+1:i*size)=(1-cos(2*pi*x/L));
    end
    U0=U0';
end

