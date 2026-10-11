function [gamma,C,info] = theorem11_bound(a,b,xi,h,L,alpha,N,epsilon)
%THEOREM11_BOUND Construct an admissible decay bound from Theorem 1.1.

a=a(:); b=b(:); xi=xi(:); h=h(:);
if any([numel(a),numel(b),numel(xi),numel(h)] ~= N)
    error('a,b,xi,h must contain N entries.');
end
if any(b<=0) || any(h<=0)
    error('This routine requires strictly positive b_j and h_j.');
end
if alpha <= N/2
    error('Theorem 1.1 requires alpha>N/2.');
end
if L >= sqrt(3)*pi/2
    error('Theorem 1.1 requires L<sqrt(3)*pi/2.');
end

epsilon_max=3*(3*pi^2-4*L^2)/(2*pi^2*L^(3/2));
if epsilon<=0 || epsilon>=epsilon_max
    error('epsilon must satisfy 0<epsilon<%g.',epsilon_max);
end

r1=(2*a-b-xi)./(L*b);
r2=(xi-b)./(L*b);
mu1_max=min([1;(2*alpha-N)/N;r1;r2]);
if mu1_max<=0
    error('No positive admissible mu1 for these parameters.');
end
mu1=0.99*mu1_max;

mu2_max=min(2*a-b-xi-mu1*L*b);
if mu2_max<=0
    error('No positive admissible mu2 for these parameters.');
end
mu2=0.99*mu2_max;

gamma_space=(3*mu1*pi^2-(2/3)*L^(3/2)*epsilon*mu1*pi^2 ...
    -4*mu1*L^2)/(8*L^2*(1+L*mu1));
gamma_delay=min(mu2./(2*h.*(mu2+xi)));
gamma=min(gamma_space,gamma_delay);
if gamma<=0
    error('The selected parameters do not produce a positive gamma.');
end

b0=min(b);
C=1+max(L*mu1,mu2/b0);

info=struct('epsilon_max',epsilon_max,'mu1_max',mu1_max,'mu1',mu1, ...
    'mu2_max',mu2_max,'mu2',mu2,'gamma_space',gamma_space, ...
    'gamma_delay',gamma_delay);
end
