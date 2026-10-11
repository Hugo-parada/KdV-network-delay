function Ufull = recover_network_state(U,N,Nx)
%RECOVER_NETWORK_STATE Recover u_{j,0},...,u_{j,Nx} for every branch.
%
% The last two values are the imposed right-boundary values
% u_{j,Nx-1}=u_{j,Nx}=0.

m = Nx-1;
Ufull = zeros(Nx+1,N);
for j = 1:N
    idx = (j-1)*m+(1:m);
    Ufull(1:m,j) = U(idx);
end
end
