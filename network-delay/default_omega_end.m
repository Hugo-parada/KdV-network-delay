function omega_end = default_omega_end(N,L)
%DEFAULT_OMEGA_END Supports used in the numerical experiments.
if N == 4
    omega_end = [L/2,L/4,L/2,L/4];
else
    omega_end = (L/2)*ones(1,N);
end
end
