# KdV equation on a star-shaped network with internal delay

MATLAB codes accompanying the numerical section of:

> H. Parada, E. Crépeau, and C. Prieur, **Delayed stabilization of the Korteweg-de Vries equation on a star-shaped network**, *Mathematics of Control, Signals, and Systems* **34** (2022), 559-605.  
> DOI: https://doi.org/10.1007/s00498-022-00319-0  
> HAL: https://hal.science/hal-03177446


## Repository layout

- `single-interval/`: nonlinear KdV on one interval.
- `network/`: nonlinear KdV on a star network without delay.
- `network-delay/`: delayed network and the scripts used for the numerical experiments.
- `common/`: routines used by the non-delayed baseline codes.
- 
Inside `network-delay/`:

- `MatDELAY.m`: matrix of the implicit KdV part, Kirchhoff condition, branchwise forward derivative and delay transport matrix.
- `omega.m`: feedback regions and coefficients.
- `initial_dataDELAY.m`: initial KdV state and initial history.
- `nonlinearpart.m`: explicit nonlinear term.
- `fullsol.m`: reconstruction of the two zero values at the exterior endpoint.
- `plotkdvdnet.m`: main simulation routine.
- `kdvnetworkDELAY.m`: direct example.
- `tests.m`-`tests4.m`: numerical experiments corresponding to Figures 2-5.

## Model represented by `network-delay`

For branch `j=1,...,N`,

```text
u_j,t + u_j,x + u_j,xxx + u_j u_j,x
      + a_j(x)u_j(t,x) + b_j(x)u_j(t-h_j,x)=0.
```

At the common node,

```text
u_j(t,0)=u_k(t,0),
sum_j u_j,xx(t,0)=-alpha*u_1(t,0)-(N/3)*u_1(t,0)^2.
```

At the exterior endpoint,

```text
u_j(t,L)=u_j,x(t,L)=0.
```

The numerical experiments use `N=4`, `L=2`, `alpha=3` and

```text
omega_1=(0,L/2),   omega_2=(0,L/4),
omega_3=(0,L/2),   omega_4=(0,L/4).
```

## Spatial discretization

Let `dx=L/Nx`. On each branch the stored unknowns are

```text
u_{j,0},...,u_{j,Nx-2}.
```

The two exterior values

```text
u_{j,Nx-1}=u_{j,Nx}=0
```

are eliminated from the unknown vector.

For interior grid points the linear KdV operator is discretized by

```text
D_x u_i=(u_{i+1}-u_{i-1})/(2*dx),
D_x^+D_x^+D_x^-u_i
=(u_{i+2}-3u_{i+1}+3u_i-u_{i-1})/dx^3.
```

The linear part is implicit, while

```text
u_i^n D_x^+u_i^n
```

is explicit.

## Common node

The first matrix row discretizes the global Kirchhoff condition directly:

```text
sum_j (u_{j,0}-2u_{j,1}+u_{j,2})/dx^2
+alpha*u_0=-(N/3)*(u_0^n)^2.
```

The remaining node rows impose

```text
u_{j,0}=u_{1,0},   j=2,...,N.
```

Consequently the common-node linear coefficient is `N/dx^2+alpha` after continuity is used.

## Delay variable

The delayed state is written as

```text
z_j(t,rho,x)=u_j(t-h_j*rho,x),
```

so that

```text
h_j z_j,t+z_j,rho=0,
z_j(t,0,x)=u_j(t,x).
```

With `drho=1/Nrho`, `Dz` in `MatDELAY.m` is the backward/upwind difference in `rho`. At each time step `plotkdvdnet.m` first computes the new memory state and then inserts its value at `rho=1` into the KdV equation. This is the same first-order splitting strategy as in the recovered code.

## Initial data

`initial_dataDELAY.m` uses

```text
u_j^0(x)=theta*(1-cos(2*pi*x/L)),
z_j^0(rho,x)=theta*(1-cos(2*pi*x/L))*cos(2*pi*rho*h_j).
```

The correction is that these functions are now evaluated on the actual grids

```text
dx=L/Nx,
drho=1/Nrho,
```

rather than on the incompatible grids `L/(Nx-2)` and `1/(Nrho-1)` used in the recovered file.

## Energy

The computed energy is the discrete counterpart of

```text
E(t)=sum_j int_0^L |u_j|^2 dx
 +sum_j h_j int_{omega_j} int_0^1 xi_j(x)|z_j|^2 d rho dx.
```

The delay term is now restricted to `omega_j`.  The vector `E` and the time vector `t` contain `Nt+1` values corresponding to

```text
t_n=n*T/Nt,   n=0,...,Nt.
```

## Numerical experiments

`tests.m` reproduces the family of feedback configurations used for Figure 2. The delay vector is

```matlab
h=[1 0.5 1 1];
```

`tests2.m` studies the delayed-feedback amplitudes/delays. 

`tests3.m` compares several amplitudes of the initial datum.

`tests4.m` compares the numerical decay with an admissible version of the bound in Theorem 1.1. 


## Running

From `network-delay/`, the simplest call is

```matlab
[E,t]=plotkdvdnet(a,b,xi,h,T,Nt,big);
```

For example,

```matlab
a=[2 4 6 8];
b=[0.5 1.5 2.5 3.5];
xi=[1 2 3 4];
h=[0.1 0.2 0.3 0.4];
[E,t]=plotkdvdnet(a,b,xi,h,10,500,0.5);
plot(t,log(E))
```
