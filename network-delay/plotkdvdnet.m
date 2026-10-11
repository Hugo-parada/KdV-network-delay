function [E,t,out] = plotkdvdnet(a,b,xi,h,T,Nt,big)
%PLOTKDVDNET Compatibility wrapper for the original numerical scripts.
[t,E,out] = solve_kdv_network_delay(a,b,xi,h,T,Nt,big);
end
