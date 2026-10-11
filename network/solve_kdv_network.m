function [t,E,out] = solve_kdv_network(N,L,alpha,T,Nt,Nx,theta)
%SOLVE_KDV_NETWORK Nonlinear KdV on a star-shaped network without delay.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));

[t,E,out] = solve_network_nodelay_core(N,L,alpha,T,Nt,Nx,theta);
end
