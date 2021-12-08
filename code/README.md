# Instructions

This directory contains code and instructions to reproduce the numerical
experiments reported in the article

> Efficient implementation of modern entropy-based discontinuous Galerkin
> methods for conservation laws

The material is structured as follows.

- The directory [`ec_performance_fluxo`](ec_performance_fluxo)
  contains all FLUXO code used in
  Section 4.1 "Baseline performance results on Cartesian and curved meshes",
  Section 6 "Gauss collocation methods and entropy projections", and
  Section 7 "More invasive optimizations".
- The directory [`Cartesian_vs_curved`](Cartesian_vs_curved)
  contains all Trixi.jl code used in
  Section 4.1 "Baseline performance results on Cartesian and curved meshes".
- The directory [`numerical_fluxes`](numerical_fluxes)
  contains all Trixi.jl code used in
  Section 4.2 "Different versions of numerical fluxes".
- The directory [`overintegration`](overintegration)
  contains all Trixi.jl code used in
  Section 5 "Comparison to overintegration".
- The directory [`Gauss`](Gauss)
  contains all Trixi.jl code used in
  Section 6 "Gauss collocation methods and entropy projections".
- The directory [`primitive_variables`](primitive_variables)
  contains all Trixi.jl code used in
  Section 7.1 "Precomputing primitive variables for the compressible Euler equations".
- The directory [`precomputed_logs`](precomputed_logs)
  contains all Trixi.jl code used in
  Section 7.2 "Precomputing primitive variables for the compressible Euler equations".
- The directory [`inlining`](inlining)
  contains all Trixi.jl code used in
  Section 7.3 "Defining options at compile time".
- The directory [`SIMD`](SIMD)
  contains all Trixi.jl code used in
  Section 8 "Explicit SIMD optimizations".

The results were obtained using Julia v1.7.0 for Trixi.jl
and the Intel Fortran Compiler 19.1.3 for FLUXO.
