function [t,E,out] = solve_kdv_network_delay(a,b,xi,h,T,Nt,theta,opts)
%SOLVE_KDV_NETWORK_DELAY Nonlinear KdV network with internal constant delays.
%
% Linear KdV and instantaneous damping are implicit; u D_x^+u is explicit;
% h_j z_t+z_rho=0 is backward Euler/upwind.

if nargin < 8 || isempty(opts)
    opts = struct();
end

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));

N = numel(a);
if N < 1 || any([numel(b),numel(xi),numel(h)] ~= N)
    error('a, b, xi and h must contain the same number of entries.');
end

L = get_option(opts,'L',2);
alpha = get_option(opts,'alpha',3);
Nx = get_option(opts,'Nx',100);
Nrho = get_option(opts,'Nrho',100);
omega_end = get_option(opts,'omega_end',default_omega_end(N,L));

if L <= 0 || T <= 0 || Nt < 1 || Nx < 4 || Nrho < 1
    error('Invalid geometry or discretization parameters.');
end
if any(h<0)
    error('Delays h_j must be nonnegative.');
end
if alpha < 0
    error('alpha must be nonnegative.');
end

m = Nx-1;
dt = T/Nt;
[a_profile,b_profile,xi_profile,xstored] = feedback_profiles( ...
    N,Nx,L,a,b,xi,omega_end);

A = build_network_matrix(N,Nx,L,dt,alpha,a_profile);
[U,Z] = initial_delay_state(N,Nx,Nrho,L,h,theta);

t = linspace(0,T,Nt+1);
E = zeros(1,Nt+1);
E(1) = network_delay_energy(U,Z,h,xi_profile,N,Nx,Nrho,L);

store_solution = (nargout >= 3);
if store_solution
    Uhist = zeros(Nx+1,N,Nt+1);
    Uhist(:,:,1) = recover_network_state(U,N,Nx);
end

for n = 1:Nt
    [Znew,Zend] = update_delay_state(U,Z,h,dt,Nrho,N,Nx);
    F = network_nonlinearity(U,N,Nx,L);

    rhs = zeros(N*m,1);
    rhs(1) = -F(1);

    for j = 1:N
        idx = (j-1)*m+(1:m);
        q = idx(2:end);
        rhs(q) = U(q)-dt*F(q)-dt*b_profile(j,2:end)'.*Zend(q);
    end

    U = A\rhs;
    Z = Znew;

    E(n+1) = network_delay_energy(U,Z,h,xi_profile,N,Nx,Nrho,L);
    if ~isfinite(E(n+1)) || any(~isfinite(U))
        warning('Non-finite values reached at t=%g.',t(n+1));
        E(n+1:end) = NaN;
        if store_solution
            Uhist(:,:,n+1:end) = NaN;
        end
        break;
    end

    if store_solution
        Uhist(:,:,n+1) = recover_network_state(U,N,Nx);
    end
end

if store_solution
    out = struct();
    out.x = linspace(0,L,Nx+1)';
    out.xstored = xstored;
    out.rho = (1:Nrho)/Nrho;
    out.U = Uhist;
    out.Z = Z;
    out.a_profile = a_profile;
    out.b_profile = b_profile;
    out.xi_profile = xi_profile;
    out.omega_end = omega_end;
    out.alpha = alpha;
else
    out = struct();
end
end

function value = get_option(opts,name,default_value)
if isfield(opts,name)
    value = opts.(name);
else
    value = default_value;
end
end
