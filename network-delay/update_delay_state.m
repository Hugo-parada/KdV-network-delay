function [Znew,Zend] = update_delay_state(U,Z,h,dt,Nrho,N,Nx)
%UPDATE_DELAY_STATE Backward Euler/upwind step for h_j z_t+z_rho=0.
%
% h_j(z_k^{n+1}-z_k^n)/dt
%   +(z_k^{n+1}-z_{k-1}^{n+1})/drho = 0,
% with inflow z_0^{n+1}=u_j^n.
%
% This keeps the first-order splitting used in the recovered code: the
% delay state is advanced first and z(rho=1) then enters the KdV step.

m = Nx-1;
drho = 1/Nrho;
q = dt/drho;
h = h(:)';

Znew = cell(1,N);
Zend = zeros(N*m,1);

for j = 1:N
    idx = (j-1)*m+(1:m);
    uj = U(idx);

    T = (h(j)+q)*speye(Nrho);
    if Nrho > 1
        T = T+sparse(2:Nrho,1:Nrho-1,-q*ones(Nrho-1,1),Nrho,Nrho);
    end

    rhs = h(j)*Z{j}';
    rhs(1,:) = rhs(1,:)+q*uj';

    Zj = (T\rhs)';
    Znew{j} = Zj;
    Zend(idx) = Zj(:,end);
end
end
