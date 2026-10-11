function [U,Z] = initial_delay_state(N,Nx,Nrho,L,h,theta)
%INITIAL_DELAY_STATE Initial state and history used in Section 5.
%
% u_j^0(x)=theta*(1-cos(2*pi*x/L)),
% z_j^0(rho,x)=u_j^0(x)*cos(2*pi*rho*h_j).

U = initial_network_state(N,Nx,L,theta);
m = Nx-1;
rho = (1:Nrho)/Nrho;
h = h(:)';

Z = cell(1,N);
for j = 1:N
    idx = (j-1)*m+(1:m);
    uj = U(idx);
    Z{j} = uj*cos(2*pi*h(j)*rho);
end
end
