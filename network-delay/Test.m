% short test of the corrected delayed-network code
clear;clc;

a=[1 1 1 1];
b=[0.5 0.5 0.5 0.5];
xi=[1 1 1 1];
h=[0.1 0.2 0.3 0.4];

[E,t]=plotkdvdnet(a,b,xi,h,0.5,20,0.2);
assert(length(E)==length(t));
assert(all(isfinite(E)));
assert(all(E>0));
disp('Delayed-network test passed.')
