function[A,Dplus]=Mat1(Nx,L,dt,alpha)
N=Nx-1;
dx=L/Nx; 


%beta1=1;
%beta2=2;
%beta3=3;
%beta4=4;
%beta5=5;
%beta6=6;
%beta7=7;

beta1=((-1/dx^3)-1/(2*dx));
beta2=(1/dt+3/dx^3);
beta3=(-3/dx^3+1/(2*dx));
beta4=1/dx^3;
beta5=(1/dx^2+alpha);
beta6=(-2/dx^2);
beta7=(1/dx^2);

A = diag(beta2*ones(1,N)) + diag(beta3*ones(1,N-1),1) + diag(beta1*ones(1,N-1),-1)+diag(beta4*ones(1,N-2),2);
A(1,1)=beta5;
A(1,2)=beta6;
A(1,3)=beta7;

Dplus=1/dx*( diag(-1*ones(1,N)) + diag(1*ones(1,N-1),1) ) ;
end

