function[U0]=inital_data1(Nx,L)
   dx2=L/(Nx-2);
    x=0:dx2:L;
    %for i=1:Nx-1;
      %  U0(i)=5;
    %end
    U0=(1-cos(2*pi*x/L));
    U0=U0';
end

