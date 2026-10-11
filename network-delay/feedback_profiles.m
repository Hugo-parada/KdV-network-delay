function [a_profile,b_profile,xi_profile,x] = feedback_profiles( ...
    N,Nx,L,a,b,xi,omega_end)
%FEEDBACK_PROFILES Piecewise-constant localized feedbacks.
%
% For the paper N=4 and
% omega_1=(0,L/2), omega_2=(0,L/4),
% omega_3=(0,L/2), omega_4=(0,L/4).

m = Nx-1;
dx = L/Nx;
x = (0:m-1)'*dx;

a = a(:)';
b = b(:)';
xi = xi(:)';
omega_end = omega_end(:)';

if any([numel(a),numel(b),numel(xi),numel(omega_end)] ~= N)
    error('a, b, xi and omega_end must contain N entries.');
end
if any(a<0) || any(b<0) || any(xi<0)
    error('Feedback and energy weights must be nonnegative.');
end

a_profile = zeros(N,m);
b_profile = zeros(N,m);
xi_profile = zeros(N,m);

for j = 1:N
    mask = (x>0) & (x<=omega_end(j)+10*eps(max(1,L)));
    a_profile(j,mask) = a(j);
    b_profile(j,mask) = b(j);
    xi_profile(j,mask) = xi(j);
end
end
