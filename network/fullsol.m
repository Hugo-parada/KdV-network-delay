function[fullU0]=fullsol(U0,N,Nx);
size=Nx-1;
fullU0=0;
for i=1:N;
    fullU0=[fullU0 ; U0((i-1)*size+1:i*size); 0 ;0];
end
fullU0=fullU0(2:end);
end