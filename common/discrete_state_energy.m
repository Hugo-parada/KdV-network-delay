function E = discrete_state_energy(U,N,Nx,L)
%DISCRETE_STATE_ENERGY Rectangle-rule approximation of sum_j int u_j^2 dx.

dx = L/Nx;
m = Nx-1;
if numel(U) ~= N*m
    error('U has incompatible size.');
end
E = dx*sum(U(:).^2);
end
