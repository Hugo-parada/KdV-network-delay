function[Aj,Bj,O]=omega(N,Nx,a,b)

size=Nx-1;
O1=zeros(1,size);
O2=zeros(1,size);
O3=zeros(1,size);
O4=zeros(1,size);

% stored index i corresponds to x=(i-1)L/Nx.
% The feedback sets are (0,L/2),(0,L/4),(0,L/2),(0,L/4).
for i=2:size;
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
end
