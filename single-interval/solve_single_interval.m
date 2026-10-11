function [t,E,out] = solve_single_interval(L,alpha,T,Nt,Nx,theta)
%SOLVE_SINGLE_INTERVAL Nonlinear KdV baseline on one interval.
%
% This is the N=1 instance of the same discretization used for the network.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'common'));

[t,E,out] = solve_network_nodelay_core(1,L,alpha,T,Nt,Nx,theta);
end
