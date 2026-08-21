# MATLAB-ODE-CPP

A C++ port of MATLAB's ODE solver suite, built for control systems research and embedded deployments.

## Declaration

This project is an independent open-source implementation of MATLAB-style ODE solver interfaces in C++.

MATLAB is a registered trademark of MathWorks. This project is not affiliated with or endorsed by MathWorks.


## Motivation

MATLAB's ODE solver suite, in particular ode45 is the default choice for simulating continuous-time dynamics during controller design. But when that controller moves from simulation to a C++ Software-In-Loop (SIL) or embedded target, the integration scheme often changes too - typically to a fixed-step Euler integrator for simplicity.

The substitution is rarely benign. Differences in step size, error control, and integration order can produce meaningfully different trajectories for the same system, particularly for fast or marginally stable dynamics. The result is a simulation-to-deployment mismatch that can silently invalidate SIL/HIL validation.

This library provides faithful C++ equivalents of MATLAB's ODE suite, so the same integration behavior carries through from MATLAB simulation to C++ deployment.


## Solvers

The following solvers have been planned to be included in this library:
- ode45: Dormand-Prince 4th/5th Order Pair
- ode23: Bogacki-Shampine 2nd/3rd Order Pair
- ode78: Runge-Kutta 7th/8th Order
- ode89: Runge-Kutta 8th/9th Order
- ode113: Variable Order Adams-Bashford-Moulton PECE Method
- ode15s: Variable 1st-5th Order Backward Differentiation Formula
- ode23s: Modified Rosenbrock Method Second Order
- ode23t: Trapezoidal Integration with Free Interpolant
- ode23tb: Trapezoidal Integration with Backward Differentiation Formula


## Status

- [x] ode45
- [x] ode23
- [x] ode78
- [x] ode89
- [x] ode113
- [x] ode15s
- [x] ode23s
- [x] ode23t
- [x] ode23tb

## License

MIT - free to use with attribution.

## Note
The code repository is part of a professional academic portfolio. Unauthorized use of code for training machine learning models or LLMs is strictly prohibited.
