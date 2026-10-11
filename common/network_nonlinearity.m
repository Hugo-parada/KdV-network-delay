function F = network_nonlinearity(U,N,Nx,L)
%NETWORK_NONLINEARITY Explicit nonlinear terms.
%
% Interior rows contain u_{j,i} D_x^+ u_{j,i}, computed branch by branch,
% so a forward difference never jumps from one edge to the next.
% The first row contains the nonlinear Kirchhoff term (N/3)u_0^2.

m = Nx-1;
dx = L/Nx;
if numel(U) ~= N*m
    error('U has incompatible size.');
end
U = U(:);
F = zeros(size(U));

for j = 1:N
    idx = (j-1)*m+(1:m);
    uj = U(idx);
    uj_plus = [uj(2:end);0]; % u_{j,Nx-1}=0
    fj = uj.*(uj_plus-uj)/dx;
    F(idx(2:end)) = fj(2:end);
end

F(1) = (N/3)*U(1)^2;
end
