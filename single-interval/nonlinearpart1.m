function[NL]=nonlinearpart1(U0,Dplus);
NL=U0.*(Dplus*U0);
NL(1)=(-1/3)*(U0(1)^2);

end