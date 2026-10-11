function A = build_network_matrix(N,Nx,L,dt,alpha,a_profile)
%BUILD_NETWORK_MATRIX Semi-implicit matrix for the nonlinear KdV network.
%
% On each branch, the stored unknowns are
%   u_{j,0},...,u_{j,Nx-2}.
% The two right boundary values u_{j,Nx-1}=u_{j,Nx}=0 are eliminated.
%
% Row 1 imposes the discrete Kirchhoff condition
%   sum_j D_xx^+ u_j(0) + alpha*u_0
%       = -(N/3)*(u_0^n)^2,
% while rows (j-1)*(Nx-1)+1, j=2,...,N impose continuity at the node.
% All other rows are backward Euler for
%   u_t + u_x + u_xxx + a(x)u
% with D_x + D_x^+ D_x^+ D_x^-.

m = Nx-1;
M = N*m;
dx = L/Nx;

if ~isequal(size(a_profile),[N,m])
    error('a_profile must have size N-by-(Nx-1).');
end
if Nx < 4
    error('Nx must be at least 4.');
end

A = spalloc(M,M,5*M);

% Kirchhoff row. Keep the N copies of u_{j,0}; the continuity rows make
% them equal, so this is exactly (N/dx^2+alpha)u_0 after elimination.
for j = 1:N
    base = (j-1)*m;
    A(1,base+1) = A(1,base+1) + 1/dx^2;
    A(1,base+2) = A(1,base+2) - 2/dx^2;
    A(1,base+3) = A(1,base+3) + 1/dx^2;
end
A(1,1) = A(1,1) + alpha;

% Continuity u_{j,0}=u_{1,0}, j=2,...,N.
for j = 2:N
    row = (j-1)*m+1;
    A(row,row) = 1;
    A(row,1) = -1;
end

c_minus = -1/dx^3 - 1/(2*dx);
c_zero  =  3/dx^3;
c_plus  = -3/dx^3 + 1/(2*dx);
c_plus2 =  1/dx^3;

% Interior evolution rows i=1,...,Nx-2 on each branch.
for j = 1:N
    base = (j-1)*m;
    for q = 2:m
        row = base+q;
        A(row,row) = 1 + dt*(c_zero+a_profile(j,q));
        A(row,base+q-1) = dt*c_minus;

        if q+1 <= m
            A(row,base+q+1) = dt*c_plus;
        end
        % If q+1=m+1, the missing value is u_{j,Nx-1}=0.

        if q+2 <= m
            A(row,base+q+2) = dt*c_plus2;
        end
        % Missing q+2 values are right-boundary zeros.
    end
end
end
