% Short consistency test for the delayed-network solver.
clear; clc;

opts = struct('Nx',30,'Nrho',20,'L',2,'alpha',3);
a = [1,1,1,1];
b = [0.5,0.5,0.5,0.5];
xi = [1,1,1,1];
h = [0.1,0.2,0.3,0.4];

[t,E,out] = solve_kdv_network_delay(a,b,xi,h,0.5,20,0.2,opts);

assert(all(isfinite(E)) && all(E>0));
node_values = squeeze(out.U(1,:,end));
assert(max(abs(node_values-node_values(1))) < 1e-9);
assert(max(abs(squeeze(out.U(end,:,end)))) < 1e-12);
assert(max(abs(squeeze(out.U(end-1,:,end)))) < 1e-12);
disp('Delayed-network test passed.');
