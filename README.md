# High-Performance CFD Development with OpenFOAM

**Team:** Ettore Cirillo, Mattia Gotti, Giulio Martella, Michele Milani, Stefano Pedretti, Daniele Piano, Federico Pinto

> ⚠️ **Repository Note:** This repository contains only the test folder and our custom user-side code. It does not include the full OpenFOAM source code, as our development methodology deliberately avoided modifying the core engine to ensure long-term maintainability.

## 📌 Overview
This project was developed for the second Module of the High Performance Scientific Computing in Aerospace Engineering course (A.Y. 2025/2026) and extends the OpenFOAM computational fluid dynamics framework to simulate high-enthalpy, hypersonic flows. Standard CFD solvers assume local thermal equilibrium, which fails in hypersonic regimes where internal energy modes (vibrational and electronic) relax much slower than translational and rotational modes. 

To solve this, we implemented a **two-temperature thermodynamic model** coupled with the **Mutation++** library. This approach allows us to accurately capture thermal non-equilibrium and finite-rate chemistry without hardcoding complex physical data into the CFD solver.

## 🛠️ Technologies
* **Frameworks & Libraries:** OpenFOAM, Mutation++
* **Languages & Parallelism:** C++, OpenMP
* **Core Concepts:** Hypersonic aerodynamics, two-temperature models, finite-rate chemistry, stiff ODE integration, and shared-memory parallelization.

## 🚀 Methodology & Key Features

* **Thermochemical Coupling:** Designed a clean API interface to query the Mutation++ backend for state-dependent properties, chemical production rates, and energy relaxation terms, effectively decoupling macroscopic transport from microscopic physics.
* **Zero-Dimensional Validation:** Developed a 0D heat bath solver to isolate and validate the thermochemical source terms against reference literature, ensuring accurate energy transfer mechanisms before multi-dimensional integration.
* **Numerical Optimization:** Implemented an adaptive time-stepping algorithm to handle the extreme numerical stiffness of the initial high-temperature transients. Added low-level C++ optimizations, such as pre-computing invariants and minimizing division operations, to reduce CPU cycles.
* **Parallel Computing:** Applied OpenMP for shared-memory parallelism to accelerate the optimized sequential kernel. Addressed thread-safety constraints within Mutation++ by employing a thread-local ownership model, successfully performing strong and weak scaling tests.
* **Non-intrusive Architecture:** Built the extension entirely in the user-space. The custom thermodynamic classes are compiled as shared object libraries and dynamically linked at runtime, preserving the integrity of the core OpenFOAM installation.
