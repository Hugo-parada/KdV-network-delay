clear all
clc


T=2;
Nt=100;  

h=[0.1 0.2 0.3 0.4];  %delay for each branch


%caso 1

a_1=[ 2 4 6 8];  %feedback for each branch
b_1=[ 0.5 1.5 2.5 3.5];  %delayed feedback for each branch 
xi_1=[1 2 3 4];


L=2;
eps=(3*(3*pi^2-4*L^2)/(2*pi^2*L^(3/2)));


mu_1=2;
for j=1:4;
   x=(2*a_1(j)-b_1(j)-xi_1(j))/(2*b_1(j));
   y=(xi_1(j)-b_1(j))/(2*b_1(j));
   z=min([1,0.5,x,y]);
   if mu_1>z;
       mu_1=z;
   end
end

mu_2=10;
for j=1:4;
   x=2*a_1(j)-b_1(j)-xi_1(j)-mu_1*2*b_1(j);
   if mu_2>x;
       mu_2=x;
   end
end


otro=1000;

for j=1:4;
   x=mu_2/(2*h(j)*(mu_2+xi_1(j)));
   if otro>x;
       otro=x;
   end
end

gamma=min((3*mu_1*pi^2-(2/3)*L^(3/2)*eps*mu_1*pi^2-mu_1*4*L^2)/(8*L^2*(1+L*mu_1)),otro);

C=(1+max(2*mu_1,mu_2*b_1(1)));


[E_1,t_1]=plotkdvdnet(a_1,b_1,xi_1,h,T,Nt,eps);

teo=log(E_1(1)*C)+(-2*gamma)*t_1;

plot(t_1,log(E_1),t_1,teo)
legend('exp','teo')



xlabel('time ') 
ylabel('ln(E)')
%set(gca, 'YScale', 'log')





