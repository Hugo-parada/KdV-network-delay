# KdV equation on a star-shaped network with internal delay

Corrected MATLAB codes accompanying the numerical section of:

> H. Parada, E. Crépeau, and C. Prieur, **Delayed stabilization of the Korteweg-de Vries equation on a star-shaped network**, *Mathematics of Control, Signals, and Systems* **34** (2022), 559-605.  
> DOI: https://doi.org/10.1007/s00498-022-00319-0  
> HAL: https://hal.science/hal-03177446

These files were reconstructed from the original simulation archive and then checked against the continuous model and the finite-difference scheme described in Section 5 of the paper.

## Repository layout

- `single-interval/`: one-interval nonlinear KdV baseline.
- `network/`: nonlinear KdV equation on a star-shaped network without delay.
- `network-delay/`: the delayed network, the experiments corresponding to Figures 2-5, and the energy/theoretical-bound routines.
- `common/`: spatial network assembly, nonlinear term, initial data and state reconstruction shared by all solvers.
- `smoke_test.m`: short consistency tests.

The PNG files already present in `network-delay/` are the figures recovered from the old archive. They were generated before the corrections listed below. Rerun `pruebas.m`-`pruebas4.m` to obtain figures from the corrected implementation.

## Model represented by the code

For branch `j=1,...,N` the controlled equation is

```text
u_j,t + u_j,x + u_j,xxx + u_j u_j,x
      + a_j(x) u_j(t,x) + b_j(x) u_j(t-h_j,x) = 0.
```

At the common vertex,

```text
u_j(t,0) = u_k(t,0),

sum_j u_j,xx(t,0)
    = -alpha*u_1(t,0) - (N/3)*u_1(t,0)^2.
```

At the exterior endpoint,

```text
u_j(t,L) = u_j,x(t,L) = 0.
```

For the numerical experiments `N=4`, `L=2`, `alpha=3` and

```text
omega_1 = (0,L/2),   omega_2 = (0,L/4),
omega_3 = (0,L/2),   omega_4 = (0,L/4).
```

## Spatial grid and network vertex

Let `dx=L/Nx`. On each branch the unknown vector stores

```text
u_{j,0},...,u_{j,Nx-2},
```

with `u_{j,Nx-1}=u_{j,Nx}=0` eliminated.

The first row is the global Kirchhoff condition

```text
sum_j (u_{j,0}-2*u_{j,1}+u_{j,2})/dx^2
    + alpha*u_0 = -(N/3)*(u_0^n)^2.
```

The next `N-1` vertex rows impose `u_{j,0}=u_{1,0}`. Hence the common-node coefficient is `N/dx^2+alpha`.

**Node-discretization note.** Section 5 of the published paper writes the forward second-derivative condition branch by branch, with coefficient `1/dx^2+alpha` and right-hand side `-(1/3)u_0^2`. Taken literally and summed over branches, that does not reproduce the global Kirchhoff condition because it multiplies `alpha*u_0` by `N`. The recovered MATLAB inherited this mismatch. The corrected code discretizes the global condition of the model directly rather than reproducing that branchwise line literally.

## KdV finite differences

For interior points:

```text
D_x u_i = (u_{i+1}-u_{i-1})/(2*dx),
D_x^+ u_i = (u_{i+1}-u_i)/dx,
D_x^+ D_x^+ D_x^- u_i
 = (u_{i+2}-3*u_{i+1}+3*u_i-u_{i-1})/dx^3.
```

The linear KdV operator and instantaneous damping are implicit. The nonlinear term `u_{j,i}^n D_x^+u_{j,i}^n` is explicit. `network_nonlinearity.m` computes forward differences separately on each branch, so the last point of branch `j` never reads the first point of branch `j+1`.

## Delay variable

For each branch,

```text
z_j(t,rho,x)=u_j(t-h_j*rho,x),
h_j z_j,t+z_j,rho=0,
z_j(t,0,x)=u_j(t,x).
```

With `drho=1/Nrho`, the transport step is

```text
h_j*(z_k^{n+1}-z_k^n)/dt
 +(z_k^{n+1}-z_{k-1}^{n+1})/drho=0.
```

The recovered code used first-order splitting: advance `z` using the current KdV state as inflow, then insert `z(rho=1)` in the KdV step. `update_delay_state.m` preserves this strategy while fixing its grid and indexing. For `h_j=0`, the same algebra makes the discrete delay state rho-constant and equal to the current state.

## Initial data

```text
u_j^0(x)=theta*(1-cos(2*pi*x/L)),
z_j^0(rho,x)=theta*(1-cos(2*pi*x/L))*cos(2*pi*rho*h_j).
```

The corrected code evaluates them on exactly the `dx=L/Nx` and `drho=1/Nrho` grids used during evolution.

## Energy

```text
E(t)=sum_j int_0^L |u_j|^2 dx
 +sum_j h_j int_{omega_j} int_0^1 xi_j(x)|z_j|^2 d rho dx.
```

`network_delay_energy.m` uses its rectangle-rule counterpart and restricts the memory term to `omega_j` through `xi_profile`. The solver returns `Nt+1` energies at the matching times `t_n=n*T/Nt`, `n=0,...,Nt`.

## Main functions

Delayed network:

```matlab
[t,E,out] = solve_kdv_network_delay(a,b,xi,h,T,Nt,theta,opts);
```

Defaults: `L=2`, `alpha=3`, `Nx=100`, `Nrho=100`.

Compatibility wrapper:

```matlab
[E,t,out] = plotkdvdnet(a,b,xi,h,T,Nt,theta);
```

Network without delay:

```matlab
[t,E,out] = solve_kdv_network(N,L,alpha,T,Nt,Nx,theta);
```

Single interval:

```matlab
[t,E,out] = solve_single_interval(L,alpha,T,Nt,Nx,theta);
```

## Numerical experiments

Run from `network-delay/`:

- `pruebas.m`: feedback configurations of Figure 2.
- `pruebas2.m`: delayed-feedback/delay comparison of Figure 3.
- `pruebas3.m`: normalized energies for `theta=0.5,2,5` as in Figure 4.
- `pruebas4.m`: numerical decay versus an admissible Theorem 1.1 bound.

For Figure 5, the theorem assumes `||U_0||_H<epsilon`; the old script used the upper bound for `epsilon` itself as the scalar amplitude `theta`. These are not the same quantity. The corrected script chooses a small `theta`, computes the discrete initial H norm and then chooses an admissible `epsilon` above it. `theorem11_bound.m` also uses strict interior values of `mu_1,mu_2` and the correct `C=1+max(L*mu_1,mu_2/b_0)`.

## Corrections made to the recovered code

1. `kdvnetworkDELAY.m` is executable; the old version reached `U0=Ui` with every definition of `Ui` commented out.
2. Initial data now use `L/Nx`, not `L/(Nx-2)`.
3. Initial histories now use `1/Nrho`, not `1/(Nrho-1)`.
4. Forward differences no longer connect consecutive branches.
5. The global Kirchhoff condition is assembled directly, giving `N/dx^2+alpha` after continuity.
6. The nonlinear term has the correct sign in the time step.
7. Memory energy is integrated only on `omega_j`.
8. Energies and times use the same `Nt+1` grid.
9. Figure 2 uses `h=[1,0.5,1,1]` as stated in Section 5; the recovered script had `h_3=0.5`.
10. Figure 3 uses `xi_3=25` when `b_3=25`; the recovered script accidentally retained `xi_3=20`.
11. The Theorem 1.1 constant uses `mu_2/b_0`, not `mu_2*b_0`.
12. `mu_1`, `mu_2` and `epsilon` are chosen strictly inside the admissible ranges.
13. The theorem's H-norm smallness parameter is no longer identified with the scalar amplitude `theta`.

## Validation

From the repository root run

```matlab
smoke_test
```

It checks finite positive energies, continuity at the central node and the two right-boundary values. The corrected algebra was also independently mirrored in Python during the cleanup to check dimensions and finite/positive energies for representative cases.

No MATLAB/Octave runtime was available in the environment used for this cleanup, so the final MATLAB execution should be performed locally with `smoke_test` before regenerating publication-quality figures.
