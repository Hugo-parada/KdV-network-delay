function[Aj,Bj]=omega(N,Nx,a,b);
O1=zeros(1,Nx-1);
O2=zeros(1,Nx-1);
O3=zeros(1,Nx-1);
O4=zeros(1,Nx-1);

for i=1:Nx-1;
    if i<=Nx/4;
        O2(i)=1;
        O4(i)=1;
        
    end
    if i<=Nx/2;
        O1(i)=1;
        O3(i)=1;
    end
end
O=[O1 O2 O3 O4];
Aj=[a(1)*O1 a(2)*O2 a(3)*O3 a(4)*O4];
Bj=[b(1)*O1 b(2)*O2 b(3)*O3 b(4)*O4];
%recall that our omega are (0,L/2) (0,L/4) (0,L/2) (0,L/4)