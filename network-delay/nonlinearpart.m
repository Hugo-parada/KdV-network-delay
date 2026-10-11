function[NL]=nonlinearpart(U0,Dplus,N);
NL=U0.*(Dplus*U0);
NL(1)=(N/3)*(U0(1)^2);

end