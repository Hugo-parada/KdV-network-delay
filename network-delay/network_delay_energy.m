function E = network_delay_energy(U,Z,h,xi_profile,N,Nx,Nrho,L)
%NETWORK_DELAY_ENERGY Discrete counterpart of the paper energy.
%
% E=sum_j int |u_j|^2 dx
%   +sum_j h_j int_{omega_j} int_0^1 xi_j(x)|z_j|^2 d rho dx.

m = Nx-1;
dx = L/Nx;
drho = 1/Nrho;
if numel(U) ~= N*m
    error('U has incompatible size.');
end

E = dx*sum(U(:).^2);
for j = 1:N
    if h(j) ~= 0
        w = xi_profile(j,:)';
        E = E+h(j)*dx*drho*sum(sum((w*ones(1,Nrho)).*(Z{j}.^2)));
    end
end
end
