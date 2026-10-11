function U = initial_network_state(N,Nx,L,theta)
%INITIAL_NETWORK_STATE Initial condition used in the numerical section.
%
% u_j^0(x)=theta*(1-cos(2*pi*x/L)).
% The stored points are the actual grid nodes x_i=i*dx, i=0,...,Nx-2.

m = Nx-1;
dx = L/Nx;
x = (0:m-1)'*dx;
u0 = theta*(1-cos(2*pi*x/L));

U = zeros(N*m,1);
for j = 1:N
    U((j-1)*m+(1:m)) = u0;
end
end
