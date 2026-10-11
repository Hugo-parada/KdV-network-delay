function[fullU0]=fullsol(U0,N,Nx)

size=Nx-1;
fullU0=[];
for j=1:N;
    fullU0=[fullU0;U0((j-1)*size+1:j*size);0;0];
end
end
