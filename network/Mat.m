function[AN,Dplus,vec]=Mat(N,Nx,L,dt,alpha)

size=Nx-1;
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

A = diag(beta2*ones(1,size)) + diag(beta3*ones(1,size-1),1) + diag(beta1*ones(1,size-1),-1)+diag(beta4*ones(1,size-2),2);
%A(1,1)=beta5;
%A(1,2)=beta6;
%A(1,3)=beta7;

A(1,1)=0;
A(1,2)=0;
A(1,3)=0;

Dplus=1/dx*(diag(-1*ones(1,N*size)) + diag(1*ones(1,N*size-1),1)) ;
AN=zeros(N*size);
vec=ones(N*size,1);
for i=1:N;
    AN((i-1)*size+1:i*size,(i-1)*size+1:i*size)=A;
    AN(1,(i-1)*size+2)=beta6;
    AN(1,(i-1)*size+3)=beta7;
    if i==1;
        AN(1,(i-1)*size+1)=beta5;
    else;
         AN((i-1)*size+1,(i-1)*size+1)=1;
         AN((i-1)*size+1,1)=-1;
    end
    vec((i-1)*size+1)=0;
    Dplus((i-1)*size+1,:)=zeros(1,N*size);
end

end

