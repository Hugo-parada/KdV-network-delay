% Quick baseline test for the network without delay.
clear; clc;

N = 3;
L = 2;
alpha = 3;
T = 1;
Nt = 40;
Nx = 40;
theta = 0.5;

[t,E,out] = solve_kdv_network(N,L,alpha,T,Nt,Nx,theta);

assert(all(isfinite(E)) && all(E>0));
node_values = squeeze(out.U(1,:,end));
assert(max(abs(node_values-node_values(1))) < 1e-9);
disp('Network-without-delay test passed.');
