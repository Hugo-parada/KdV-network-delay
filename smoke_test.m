% Short consistency tests for the corrected repository.
clear;clc;

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'common'));
addpath(fullfile(root,'single-interval'));
addpath(fullfile(root,'network'));
addpath(fullfile(root,'network-delay'));

[t1,E1,out1]=solve_single_interval(2,2,0.2,10,20,0.1);
assert(length(t1)==length(E1));
assert(all(isfinite(E1)) && all(E1>0));

[t2,E2,out2]=solve_kdv_network(3,2,2,0.2,10,20,0.1);
assert(length(t2)==length(E2));
assert(all(isfinite(E2)) && all(E2>0));
node_values=squeeze(out2.U(1,:,end));
assert(max(abs(node_values-node_values(1)))<1e-9);

a=[1 1 1 1];
b=[0.5 0.5 0.5 0.5];
xi=[1 1 1 1];
h=[0.1 0.2 0.3 0.4];
[E3,t3]=plotkdvdnet(a,b,xi,h,0.5,20,0.2);
assert(length(t3)==length(E3));
assert(all(isfinite(E3)) && all(E3>0));

disp('All smoke tests passed.')
