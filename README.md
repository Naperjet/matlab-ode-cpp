# C++ High-Fidelity ODE Suite & Zero-Overhead Linear Algebra Library

## Abstract
This repository contains a custom C++ Ordinary Differential Equation (ODE) solver suite and a purpose-built linear algebra library, designed to bridge the gap between high-level prototyping and deterministic Software-in-the-Loop (SIL) pipelines for aerospace Guidance, Navigation, and Control (GNC). By implementing explicit and implicit Runge-Kutta methods with strict memory and precision controls, this suite provides bit-level reproducibility and real-time performance guarantees absent in standard auto-generated MATLAB-to-C++ pipelines.

Visit GitLab for actual development timeline and standalone files

GitLab Link:[https://www.gitlab.com/naperjet/matlab-ode-cpp]

## The Core Problem
In aerospace GNC development, dynamic models are typically prototyped in MATLAB/Simulink. Transitioning these models to C++ for SIL or Hardware-in-the-Loop (HIL) integration often introduces pipeline inconsistencies, primarily stemming from dynamic memory allocation and floating-point corruption.

High-order Runge-Kutta methods rely on Butcher tableau coefficients that require 27 to 35 digits of precision to maintain stability. When naively ported to standard C++ double precision (15-17 digits), truncation of these coefficients introduces systemic numerical error, degrading solver convergence and threatening the validity of the SIL simulation. This suite mitigates that corruption through mixed-precision arithmetic and rigorous algorithmic capping.

## Technical Architecture

### High Precision Solvers
To address the truncation of Butcher tableau coefficients, the internal arithmetic is executed using long double (extended precision), while the state vector output is constrained to standard double for SIL interface compatibility. Furthermore, ode113 is strictly capped at 7th order for standard double-precision output, preventing the propagation of coefficient truncation noise inherent to higher-order variable-step methods.

- **ode78:** 7(8) explicit Runge-Kutta. Optimized for highly non-stiff, smooth dynamical systems requiring high accuracy over long integration horizons.
- **ode89:** 8(9) explicit Runge-Kutta. Reserved for ultra-high precision requirements where the computational overhead is justified by extreme error tolerances.
- **ode113:** Variable-order Adams-Bashforth-Moulton predictor-corrector. Capped at 7th order for standard double-precision environments to ensure numerical stability.
- **ode113_long:** Extended-precision implementation of the Adams-Bashforth-Moulton method, utilizing long double natively for systems requiring higher-order variable-step integration.

### Stiff and Implicit Solvers
Aerospace GNC models frequently encounter stiff dynamics (e.g., tightly coupled actuator dynamics, orbital perturbations, or localized flex modes). These solvers utilize implicit formulations requiring exact Jacobian evaluations and Newton-Raphson iterations at every integration step.

- **ode15s:** Variable-order solver based on Backward Differentiation Formulas (BDF) and Numerical Differentiation Formulas (NDF). The standard workhorse for highly stiff, large-scale dynamical systems where explicit methods fail to converge.
- **ode23s:** Rosenbrock method (2nd order). Modified for stiff systems where exact Jacobian evaluation is feasible. Optimal for moderately stiff systems with tight real-time deadlines.
- **ode23t:**  Trapezoidal rule (2nd order). Implements implicit integration without numerical damping, ensuring strict conservation of energy in marginally stiff mechanical systems.
- **ode23tb:** TR-BDF2 (2nd order). Combines the trapezoidal rule with the second-order backward differentiation formula. Provides robust L-stability, making it highly suitable for highly stiff, non-linear GNC plant models.

### The Workhorse Solvers
Most plant dynamics ode equations donot feature highly stiff or high precision requiring long double levels of precision. The workhorse solvers solve explicit, smooth and medium precision ode equations and are meant to be fast and easy to call.

- **ode45:** Dormand-Prince 4(5) explicit Runge-Kutta. The foundational solver for non-stiff, smooth GNC plant models and 6-DOF trajectory propagation.
- **ode23:** Bogacki-Shampine 3(2) explicit Runge-Kutta. The lightweight workhorse, optimized for real-time loops requiring crude tolerances and low computational overhead.


## Performance and Use Cases
This architecture is explicitly designed for computationally constrained, real-time SIL environments where jitter is unacceptable.

- **Use Case - 6-DOF Aerospace Simulation:** High-fidelity missile or spacecraft 6-DOF simulations where tightly coupled aerodynamics/flex-body dynamics require ode23tb or ode15s for stability, while the guidance loops demand deterministic memory.
- **Use Case - Bit-Exact Regression Testing:** Environments where MATLAB-generated C++ code fails to meet bit-exactness criteria against legacy Fortran baselines due to double-precision Butcher coefficient truncation. The mixed-precision approach guarantees identical trajectory propagation.
- **Performance Profile:** By eliminating mid-integration heap allocations and enforcing cache-friendly pointer arithmetic, the solvers operate with strictly bounded execution times per integration step, making them suitable for hard real-time OS (RTOS) integration.

# License and Disclaimer
This project is licensed under the MIT License.

**AI Training Disclaimer:** This codebase is part of a professional academic portfolio. Unauthorized use of this repository, its code, or its documentation for training machine learning models, large language models (LLMs), or other artificial intelligence systems is strictly prohibited without explicit written consent.
