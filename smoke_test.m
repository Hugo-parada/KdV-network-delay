%SMOKE_TEST Short consistency tests for the corrected repository.
clear; clc;

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'common'));
addpath(fullfile(root,'single-interval'));
addpath(fullfile(root,'network'));
addpath(fullfile(root,'network-delay'));

[t1,E1,out1] = solve_single_interval(2,2,0.2,10,20,0.1);
assert(numel(t1)==numel(E1));
assert(all(isfinite(E1)) && all(E1>0));
right1 = squeeze(out1.U(end,1,:));
right2 = squeeze(out1.U(end-1,1,:));
assert(max(abs(right1))<1e-12);
assert(max(abs(right2))<1e-12);

[t2,E2,out2] = solve_kdv_network(3,2,2,0.2,10,20,0.1);
assert(numel(t2)==numel(E2));
assert(all(isfinite(E2)) && all(E2>0));
nodes = squeeze(out2.U(1,:,end));
assert(max(abs(nodes-nodes(1)))<1e-9);

opts = struct('Nx',24,'Nrho',16,'L',2,'alpha',3);
a = [1,1,1,1];
b = [0.5,0.5,0.5,0.5];
xi = [1,1,1,1];
h = [0.1,0.2,0.3,0.4];
[t3,E3,out3] = solve_kdv_network_delay(a,b,xi,h,0.2,10,0.1,opts);
assert(numel(t3)==numel(E3));
assert(all(isfinite(E3)) && all(E3>0));
nodes = squeeze(out3.U(1,:,end));
assert(max(abs(nodes-nodes(1)))<1e-9);
assert(max(abs(squeeze(out3.U(end,:,end))))<1e-12);
assert(max(abs(squeeze(out3.U(end-1,:,end))))<1e-12);

disp('All smoke tests passed.');
