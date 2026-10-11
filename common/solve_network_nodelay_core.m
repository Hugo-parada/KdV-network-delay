function [t,E,out] = solve_network_nodelay_core(N,L,alpha,T,Nt,Nx,theta)
%SOLVE_NETWORK_NODELAY_CORE Nonlinear KdV on a finite star network.
%
% The linear KdV operator is implicit and u D_x^+u is explicit.

if N < 1 || Nx < 4 || Nt < 1 || L <= 0 || T <= 0
    error('Invalid discretization or geometry parameters.');
end
if alpha < 0
    error('alpha must be nonnegative.');
end

m = Nx-1;
dt = T/Nt;
a_profile = zeros(N,m);
A = build_network_matrix(N,Nx,L,dt,alpha,a_profile);
U = initial_network_state(N,Nx,L,theta);

t = linspace(0,T,Nt+1);
E = zeros(1,Nt+1);
E(1) = discrete_state_energy(U,N,Nx,L);

store_solution = (nargout >= 3);
if store_solution
    Uhist = zeros(Nx+1,N,Nt+1);
    Uhist(:,:,1) = recover_network_state(U,N,Nx);
end

for n = 1:Nt
    F = network_nonlinearity(U,N,Nx,L);
    rhs = zeros(N*m,1);

    % Nonlinear Kirchhoff condition at the new linear time level.
    rhs(1) = -F(1);

    % Continuity rows have zero right-hand side.
    for j = 1:N
        idx = (j-1)*m+(1:m);
        rhs(idx(2:end)) = U(idx(2:end))-dt*F(idx(2:end));
    end

    U = A\rhs;
    E(n+1) = discrete_state_energy(U,N,Nx,L);

    if store_solution
        Uhist(:,:,n+1) = recover_network_state(U,N,Nx);
    end
end

if store_solution
    out = struct();
    out.x = linspace(0,L,Nx+1)';
    out.U = Uhist;
    out.N = N;
    out.L = L;
    out.alpha = alpha;
else
    out = struct();
end
end
