# KdV equation on a star-shaped network with internal delay

MATLAB codes accompanying the numerical simulations in:

> H. Parada, E. Crépeau, and C. Prieur, **Delayed stabilization of the Korteweg–de Vries equation on a star-shaped network**, *Mathematics of Control, Signals, and Systems* **34** (2022), 559–605.

- DOI: [10.1007/s00498-022-00319-0](https://doi.org/10.1007/s00498-022-00319-0)
- Article: [Springer](https://link.springer.com/article/10.1007/s00498-022-00319-0)

## Contents

- `single-interval/`: baseline nonlinear KdV simulation on one bounded interval.
- `network/`: nonlinear KdV simulation on a finite star-shaped network without delay.
- `network-delay/`: star-shaped network with localized instantaneous and delayed internal feedbacks, including the scripts used for the numerical experiments and the generated figures.

## Numerical method

The spatial derivatives are approximated by finite differences. The linear KdV part is treated implicitly and the nonlinear term is evaluated explicitly. The delayed term is represented by the transport variable

```text
z_j(t,x,rho) = u_j(t-rho*h_j,x),   0 < rho < 1,
```

which is also discretized by finite differences. The code tracks the Lyapunov energy containing both the KdV state and the delay state.

## Main scripts

### Single interval

Run `single-interval/kdv1.m`.

### Network without delay

Run `network/kdvnetwork.m`.

### Network with delay

The main simulation routines are:

- `network-delay/kdvnetworkDELAY.m`
- `network-delay/Test.m`
- `network-delay/plotkdvdnet.m`

The figure scripts are:

- `network-delay/pruebas.m`: comparison of instantaneous and delayed damping coefficients.
- `network-delay/pruebas2.m`: influence of the delays and delayed feedback amplitudes.
- `network-delay/pruebas3.m`: comparison of several initial-data amplitudes.
- `network-delay/pruebas4.m`: numerical decay compared with the theoretical exponential estimate.

Run each script from its own directory so that MATLAB finds the associated functions.
