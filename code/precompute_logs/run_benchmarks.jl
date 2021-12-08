
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using DelimitedFiles, LinearAlgebra, Statistics
using BenchmarkTools
using MuladdMacro
using StaticArrays
using StrideArrays
using Trixi

import Trixi: flux_ranocha, ln_mean, inv_ln_mean

struct PrimitiveVariablesCompressibleEuler2D{T<:Real} <: FieldVector{4, T}
  rho::T
  v1 ::T
  v2 ::T
  p  ::T
end

@inline function Trixi.cons2prim(u::PrimitiveVariablesCompressibleEuler2D,
                                 equations::CompressibleEulerEquations2D)
  return u
end

struct PrimitiveVariablesCompressibleEuler3D{T<:Real} <: FieldVector{5, T}
  rho::T
  v1 ::T
  v2 ::T
  v3 ::T
  p  ::T
end

@inline function Trixi.cons2prim(u::PrimitiveVariablesCompressibleEuler3D,
                                 equations::CompressibleEulerEquations3D)
  return u
end


@inline function Trixi.get_node_vars(u::AbstractArray{<:StaticVector},
                                     equations, dg::DG, indices...)
  return u[indices...]
end


function cons2prim!(prim::AbstractVector{<:StaticVector},
                    u::AbstractArray{<:Real, 2},
                    equations, solver)
  for i in eachindex(prim)
    u_local = Trixi.get_node_vars(u, equations, solver, i)
    prim_local = cons2prim(u_local, equations)
    prim[i] = prim_local
  end
end

function cons2prim!(_prim::AbstractVector{<:PrimitiveVariablesCompressibleEuler2D},
                    u::AbstractArray{<:Real, 2},
                    equations::CompressibleEulerEquations2D, solver)
  prim = reinterpret(reshape, eltype(eltype(_prim)), _prim)

  Trixi.@turbo for i in Trixi.indices((prim, u), 2)
    rho    = u[1, i]
    rho_v1 = u[2, i]
    rho_v2 = u[3, i]
    rho_e  = u[4, i]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))

    prim[1, i] = rho
    prim[2, i] = v1
    prim[3, i] = v2
    prim[4, i] = p
  end
end

function cons2prim!(_prim::AbstractVector{<:PrimitiveVariablesCompressibleEuler3D},
                    u::AbstractArray{<:Real, 2},
                    equations::CompressibleEulerEquations3D, solver)
  prim = reinterpret(reshape, eltype(eltype(_prim)), _prim)

  Trixi.@turbo for i in Trixi.indices((prim, u), 2)
    rho    = u[1, i]
    rho_v1 = u[2, i]
    rho_v2 = u[3, i]
    rho_v3 = u[4, i]
    rho_e  = u[5, i]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    prim[1, i] = rho
    prim[2, i] = v1
    prim[3, i] = v2
    prim[4, i] = v3
    prim[5, i] = p
  end
end


struct PrecomputedVariablesCompressibleEulerEquations2D{T<:Real} <: FieldVector{6, T}
    rho::T
    v1::T
    v2::T
    p::T
    log_rho::T
    log_p::T
end

struct PrecomputedVariablesCompressibleEulerEquations3D{T<:Real} <: FieldVector{7, T}
    rho::T
    v1::T
    v2::T
    v3::T
    p::T
    log_rho::T
    log_p::T
end


@muladd begin

"""
    ln_mean(x, y)

Compute the logarithmic mean

    ln_mean(x, y) = (y - x) / (log(y) - log(x)) = (y - x) / log(y / x)

Problem: The formula above has a removable singularity at `x == y`. Thus,
some care must be taken to implement it correctly without problems or loss
of accuracy when `x ≈ y`. Here, we use the approach proposed by
Ismail and Roe (2009).
Set ξ = y / x. Then, we have

    (y - x) / log(y / x) = (x + y) / log(ξ) * (ξ - 1) / (ξ + 1)

Set f = (ξ - 1) / (ξ + 1) = (y - x) / (x + y). Then, we use the expansion

    log(ξ) = 2 * f * (1 + f^2 / 3 + f^4 / 5 + f^6 / 7) + O(ξ^9)

Inserting the first few terms of this expansion yields

    (y - x) / log(ξ) ≈ (x + y) * f / (2 * f * (1 + f^2 / 3 + f^4 / 5 + f^6 / 7))
                        = (x + y) / (2 + 2/3 * f^2 + 2/5 * f^4 + 2/7 * f^6)

Since divisions are usually more expensive on modern hardware than
multiplications (Agner Fog), we try to avoid computing two divisions. Thus,
we use

    f^2 = (y - x)^2 / (x + y)^2
        = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y)

Given ε = 1.0e-4, we use the following algorithm.

    if f^2 < ε
        # use the expansion above
    else
        # use the direct formula (y - x) / log(y / x)
    end

# References
- Ismail, Roe (2009).
    Affordable, entropy-consistent Euler flux functions II: Entropy production at shocks.
    [DOI: 10.1016/j.jcp.2009.04.021](https://doi.org/10.1016/j.jcp.2009.04.021)
- Agner Fog.
    Lists of instruction latencies, throughputs and micro-operation breakdowns
    for Intel, AMD, and VIA CPUs.
    https://www.agner.org/optimize/instruction_tables.pdf
"""
@inline function Trixi.ln_mean(x, y, log_x, log_y)
    epsilon_f2 = 1.0e-4
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return (x + y) / @evalpoly(f2, 2, 2/3, 2/5, 2/7)
    else
        return (y - x) / (log_y - log_x) # log(y / x)
    end
end

"""
    inv_ln_mean(x, y)

Compute the inverse `1 / ln_mean(x, y)` of the logarithmic mean
[`ln_mean`](@ref).

This function may be used to increase performance where the inverse of the
logarithmic mean is needed, by replacing a (slow) division by a (fast)
multiplication.
"""
@inline function Trixi.inv_ln_mean(x, y, log_x, log_y)
    epsilon_f2 = 1.0e-4
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        return @evalpoly(f2, 2, 2/3, 2/5, 2/7) / (x + y)
    else
        return (log_y - log_x) / (y - x)
    end
end


"""
    flux_ranocha(u_ll, u_rr, orientation_or_normal_direction,
                 equations::CompressibleEulerEquations2D)

Entropy conserving and kinetic energy preserving two-point flux by
- Hendrik Ranocha (2018)
  Generalised Summation-by-Parts Operators and Entropy Stability of Numerical Methods
  for Hyperbolic Balance Laws
  [PhD thesis, TU Braunschweig](https://cuvillier.de/en/shop/publications/7743)
See also
- Hendrik Ranocha (2020)
  Entropy Conserving and Kinetic Energy Preserving Numerical Methods for
  the Euler Equations Using Summation-by-Parts Operators
  [Proceedings of ICOSAHOM 2018](https://doi.org/10.1007/978-3-030-39647-3_42)
"""
@inline function flux_ranocha(u_ll::PrecomputedVariablesCompressibleEulerEquations2D,
                              u_rr::PrecomputedVariablesCompressibleEulerEquations2D,
                              orientation::Integer, equations::CompressibleEulerEquations2D)
  # Unpack left and right state
  rho_ll, v1_ll, v2_ll, p_ll, log_rho_ll, log_p_ll = u_ll
  rho_rr, v1_rr, v2_rr, p_rr, log_rho_rr, log_p_rr = u_rr

  # Compute the necessary mean values
  rho_mean = ln_mean(rho_ll, rho_rr, log_rho_ll, log_rho_rr)
  # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
  # in exact arithmetic since
  #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
  #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
  inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll, log_rho_ll + log_p_rr, log_rho_rr + log_p_ll)
  v1_avg = 0.5 * (v1_ll + v1_rr)
  v2_avg = 0.5 * (v2_ll + v2_rr)
  p_avg  = 0.5 * (p_ll + p_rr)
  velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr)

  # Calculate fluxes depending on orientation
  if orientation == 1
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = f1 * v2_avg
    f4 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)
  else
    f1 = rho_mean * v2_avg
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg + p_avg
    f4 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v2_rr + p_rr*v2_ll)
  end

  return SVector(f1, f2, f3, f4)
end


@inline function flux_ranocha(u_ll::PrecomputedVariablesCompressibleEulerEquations3D,
                              u_rr::PrecomputedVariablesCompressibleEulerEquations3D,
                              orientation::Integer, equations::CompressibleEulerEquations3D)
  # Unpack left and right state
  rho_ll, v1_ll, v2_ll, v3_ll, p_ll, log_rho_ll, log_p_ll = u_ll
  rho_rr, v1_rr, v2_rr, v3_rr, p_rr, log_rho_rr, log_p_rr = u_rr

  # Compute the necessary mean values
  rho_mean = ln_mean(rho_ll, rho_rr, log_rho_ll, log_rho_rr)
  # Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
  # in exact arithmetic since
  #     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
  #   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
  inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll, log_rho_ll + log_p_rr, log_rho_rr + log_p_ll)
  v1_avg = 0.5 * (v1_ll + v1_rr)
  v2_avg = 0.5 * (v2_ll + v2_rr)
  v3_avg = 0.5 * (v3_ll + v3_rr)
  p_avg  = 0.5 * (p_ll + p_rr)
  velocity_square_avg = 0.5 * (v1_ll*v1_rr + v2_ll*v2_rr + v3_ll*v3_rr)

  # Calculate fluxes depending on orientation
  if orientation == 1
    f1 = rho_mean * v1_avg
    f2 = f1 * v1_avg + p_avg
    f3 = f1 * v2_avg
    f4 = f1 * v3_avg
    f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)
  elseif orientation == 2
    f1 = rho_mean * v2_avg
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg + p_avg
    f4 = f1 * v3_avg
    f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v2_rr + p_rr*v2_ll)
  else # orientation == 3
    f1 = rho_mean * v3_avg
    f2 = f1 * v1_avg
    f3 = f1 * v2_avg
    f4 = f1 * v3_avg + p_avg
    f5 = f1 * ( velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one ) + 0.5 * (p_ll*v3_rr + p_rr*v3_ll)
  end

  return SVector(f1, f2, f3, f4, f5)
end

end # muladd


@inline function Trixi.get_node_vars(u::AbstractArray{<:StaticVector},
                                     equations, dg::DG, indices...)
  return u[indices...]
end


function cons2prim!(_prim::AbstractVector{<:PrecomputedVariablesCompressibleEulerEquations2D},
                    u::AbstractArray{<:Real, 2},
                    equations::CompressibleEulerEquations2D, solver)
  prim = reinterpret(reshape, eltype(eltype(_prim)), _prim)

  Trixi.@turbo for i in Trixi.indices((prim, u), 2)
    rho    = u[1, i]
    rho_v1 = u[2, i]
    rho_v2 = u[3, i]
    rho_e  = u[4, i]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))

    prim[1, i] = rho
    prim[2, i] = v1
    prim[3, i] = v2
    prim[4, i] = p
    prim[5, i] = log(rho)
    prim[6, i] = log(p)
  end
end

function cons2prim!(_prim::AbstractVector{<:PrecomputedVariablesCompressibleEulerEquations3D},
                    u::AbstractArray{<:Real, 2},
                    equations::CompressibleEulerEquations3D, solver)
  prim = reinterpret(reshape, eltype(eltype(_prim)), _prim)

  Trixi.@turbo for i in Trixi.indices((prim, u), 2)
    rho    = u[1, i]
    rho_v1 = u[2, i]
    rho_v2 = u[3, i]
    rho_v3 = u[4, i]
    rho_e  = u[5, i]

    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    v3 = rho_v3 / rho
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

    prim[1, i] = rho
    prim[2, i] = v1
    prim[3, i] = v2
    prim[4, i] = v3
    prim[5, i] = p
    prim[6, i] = log(rho)
    prim[7, i] = log(p)
  end
end


## set up benchmark code
"""
    initial_condition_isentropic_vortex(x, t, equations)

The classical isentropic vortex test case of
- Chi-Wang Shu (1997)
  Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory
  Schemes for Hyperbolic Conservation Laws.
  [NASA/CR-97-206253](https://ntrs.nasa.gov/citations/19980007543)
"""
function initial_condition_isentropic_vortex(x, t, equations::CompressibleEulerEquations2D)
  ϱ0 = 1.0               # background density
  v0 = SVector(1.0, 1.0) # background velocity
  p0 = 10.0              # background pressure
  ε = 20.0               # vortex strength
  L = 10.0               # size of the domain per coordinate direction

  T0 = p0 / ϱ0           # background temperature
  γ = equations.gamma    # ideal gas constant

  x0 = v0 * t            # current center of the vortex
  dx = vortex_center.(x - x0, L)
  r2 = sum(abs2, dx)

  # perturbed primitive variables
  T = T0 - (γ - 1) * ε^2 / (8 * γ * π^2) * exp(1 - r2)
  v = v0 + ε / (2 * π) * exp(0.5 * (1 - r2)) * SVector(-dx[2], dx[1])
  ϱ = ϱ0 * (T / T0)^(1 / (γ - 1))
  p = ϱ * T

  return prim2cons(SVector(ϱ, v..., p), equations)
end

function initial_condition_isentropic_vortex(x, t, equations::CompressibleEulerEquations3D)
  ϱ0 = 1.0               # background density
  v0 = SVector(1.0, 1.0, 0.0) # background velocity
  p0 = 10.0              # background pressure
  ε = 20.0               # vortex strength
  L = 10.0               # size of the domain per coordinate direction

  T0 = p0 / ϱ0           # background temperature
  γ = equations.gamma    # ideal gas constant

  x0 = v0 * t            # current center of the vortex
  dx = vortex_center.(x - x0, L)
  r2 = sum(abs2, dx)

  # perturbed primitive variables
  T = T0 - (γ - 1) * ε^2 / (8 * γ * π^2) * exp(1 - r2)
  v = v0 + ε / (2 * π) * exp(0.5 * (1 - r2)) * SVector(-dx[2], dx[1], 0.0)
  ϱ = ϱ0 * (T / T0)^(1 / (γ - 1))
  p = ϱ * T

  return prim2cons(SVector(ϱ, v..., p), equations)
end

vortex_center(x, L) = mod(x + L/2, L) - L/2


function split_form_prim_kernel!(du, u_prim_logs, u_cons,
                                 element, mesh,
                                 nonconservative_terms, equations,
                                 volume_flux::VF, dg, cache) where {VF}

  cons2prim!(vec(u_prim_logs), reshape(u_cons, nvariables(equations), :),
             equations, dg)

  Trixi.split_form_kernel!(du, u_prim_logs, element, mesh,
                           nonconservative_terms, equations,
                           volume_flux, dg, cache)

  return nothing
end


# 2D version
function run_benchmarks_2d(polydeg, numerical_flux, initial_condition)
  println("Benchmarking ", numerical_flux, " in 2D with polydeg=", polydeg)

  # setup basic code
  equations = CompressibleEulerEquations2D(1.4)

  coordinates_min = (-5.0, -5.0)
  coordinates_max = ( 5.0,  5.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=0, n_cells_max=100_000, periodicity=true)

  solver = DGSEM(polydeg=polydeg, surface_flux=numerical_flux,
                 volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))

  semi = SemidiscretizationHyperbolic(mesh, equations,
                                      initial_condition, solver)

  u_ode  = compute_coefficients(0.0, semi)
  du_ode = similar(u_ode)
  GC.@preserve u_ode du_ode begin
    u  = Trixi.wrap_array(u_ode,  semi)
    du = Trixi.wrap_array(du_ode, semi)

    u_prim = StrideArray(undef,
                         PrimitiveVariablesCompressibleEuler2D{eltype(u)},
                         Base.tail(StrideArrays.ArrayInterface.size(u))...)

    u_prim_logs = StrideArray(undef,
                              PrecomputedVariablesCompressibleEulerEquations2D{eltype(u)},
                              Base.tail(StrideArrays.ArrayInterface.size(u))...)
    element = 1

    # benchmark flux differencing volume terms based on conservative variables
    fill!(du, zero(eltype(du)))
    Trixi.split_form_kernel!(
      du, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_cons = copy(du)

    sleep(0.5)
    result_cons = @benchmark Trixi.split_form_kernel!(
      $du, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_cons)

    # benchmark flux differencing volume terms based on primitive variables
    fill!(du, zero(eltype(du)))
    split_form_prim_kernel!(
      du, u_prim, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_prim = copy(du)
    @assert du_cons ≈ du_prim

    sleep(0.5)
    result_prim = @benchmark split_form_prim_kernel!(
      $du, $u_prim, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_prim)

    # benchmark both primitive variables and precomputed logarithms
    fill!(du, zero(eltype(du)))
    split_form_prim_kernel!(
      du, u_prim_logs, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_prim_logs = copy(du)
    @assert du_cons ≈ du_prim_logs

    sleep(0.5)
    result_prim_logs = @benchmark split_form_prim_kernel!(
      $du, $u_prim_logs, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_prim_logs)
  end
  println()

  return result_cons, result_prim, result_prim_logs
end


# 3D version
function run_benchmarks_3d(polydeg, numerical_flux, initial_condition)
  println("Benchmarking ", numerical_flux, " in 3D with polydeg=", polydeg)

  # setup basic code
  equations = CompressibleEulerEquations3D(1.4)

  coordinates_min = (-5.0, -5.0, -5.0)
  coordinates_max = ( 5.0,  5.0,  5.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=0, n_cells_max=100_000, periodicity=true)

  solver = DGSEM(polydeg=polydeg, surface_flux=numerical_flux,
                 volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))

  semi = SemidiscretizationHyperbolic(mesh, equations,
                                      initial_condition, solver)

  u_ode  = compute_coefficients(0.0, semi)
  du_ode = similar(u_ode)
  GC.@preserve u_ode du_ode begin
    u  = Trixi.wrap_array(u_ode,  semi)
    du = Trixi.wrap_array(du_ode, semi)

    u_prim = StrideArray(undef,
                         PrimitiveVariablesCompressibleEuler3D{eltype(u)},
                         Base.tail(StrideArrays.ArrayInterface.size(u))...)

    u_prim_logs = StrideArray(undef,
                         PrecomputedVariablesCompressibleEulerEquations3D{eltype(u)},
                         Base.tail(StrideArrays.ArrayInterface.size(u))...)
    element = 1

    # benchmark flux differencing volume terms based on conservative variables
    fill!(du, zero(eltype(du)))
    Trixi.split_form_kernel!(
      du, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_cons = copy(du)

    sleep(0.5)
    result_cons = @benchmark Trixi.split_form_kernel!(
      $du, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_cons)

    # benchmark flux differencing volume terms based on primitive variables
    fill!(du, zero(eltype(du)))
    split_form_prim_kernel!(
      du, u_prim, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_prim = copy(du)
    @assert du_cons ≈ du_prim

    sleep(0.5)
    result_prim = @benchmark split_form_prim_kernel!(
      $du, $u_prim, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_prim)

    # flux differencing volume terms based on primitive variables and precomputed logarithms
    fill!(du, zero(eltype(du)))
    split_form_prim_kernel!(
      du, u_prim_logs, u,
      element, mesh,
      (Trixi.have_nonconservative_terms(equations)), equations,
      numerical_flux, solver, (semi.cache))
    du_prim_logs = copy(du)
    @assert du_cons ≈ du_prim_logs

    sleep(0.5)
    result_prim_logs = @benchmark split_form_prim_kernel!(
      $du, $u_prim_logs, $u,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver, $(semi.cache)) setup=(fill!($du, zero(eltype($du))))
    display(result_prim_logs)
  end
  println()

  return result_cons, result_prim, result_prim_logs
end

# general code
Base.ndims(::typeof(run_benchmarks_2d)) = 2
Base.ndims(::typeof(run_benchmarks_3d)) = 3

function run_benchmarks(run_benchmarks_inner, polydegs, numerical_flux, initial_condition, icname="")
  results = zeros(length(polydegs), 1)
  results .= polydegs
  header = "# Polydeg"

  # flux differencing
  println("#"^80)
  println("flux differencing\n")
  means_cons = zeros(length(polydegs))
  stds_cons  = zero(means_cons)
  means_prim = zero(means_cons)
  stds_prim  = zero(means_cons)
  means_prim_logs = zero(means_cons)
  stds_prim_logs  = zero(means_cons)
  for (i, polydeg) in enumerate(polydegs)
    result_cons, result_prim, result_prim_logs = run_benchmarks_inner(polydeg, numerical_flux, initial_condition)
    means_cons[i] = time(mean(result_cons)) * 1.0e-9
    stds_cons[i]  = time(std(result_cons))  * 1.0e-9
    means_prim[i] = time(mean(result_prim)) * 1.0e-9
    stds_prim[i]  = time(std(result_prim))  * 1.0e-9
    means_prim_logs[i] = time(mean(result_prim_logs)) * 1.0e-9
    stds_prim_logs[i]  = time(std(result_prim_logs))  * 1.0e-9
  end
  results = hcat(results, means_cons, stds_cons, means_prim, stds_prim, means_prim_logs, stds_prim_logs)
  header = header * "\tConservative-variables-mean\tstd" *
                    "\tPrimitive-variables-mean\tstd" * 
                    "\tPrecomputed-log-variables-mean\tstd"

  # write results to a file
  open(joinpath(@__DIR__, "precomputed_$(ndims(run_benchmarks_inner))D_$(numerical_flux)_$icname.dat"), "w") do io
    println(io, header)
    writedlm(io, results)
  end

  return results
end


## run benchmarks
polydegs = 3:15
run_benchmarks(run_benchmarks_2d, polydegs, flux_ranocha, initial_condition_isentropic_vortex, "vortex")
run_benchmarks(run_benchmarks_3d, polydegs, flux_ranocha, initial_condition_isentropic_vortex, "vortex")

initial_condition_random(x, t, equations::CompressibleEulerEquations2D) = prim2cons(SVector{4}(rand(), .1, .2, rand()), equations)
initial_condition_random(x, t, equations::CompressibleEulerEquations3D) = prim2cons(SVector{5}(rand(), .1, .2, .3, rand()), equations)

run_benchmarks(run_benchmarks_2d, polydegs, flux_ranocha, initial_condition_random, "random")
run_benchmarks(run_benchmarks_3d, polydegs, flux_ranocha, initial_condition_random, "random")

function initial_condition_sine(x, t, equations::CompressibleEulerEquations2D) 
    @unpack gamma = equations
    rho = 2.0 + sin(pi * x[1] / 5.0) * sin(pi * x[2] / 5.0)
    return prim2cons(SVector{4}(rho, .1, .2, rho^gamma), equations)
end
function initial_condition_sine(x, t, equations::CompressibleEulerEquations3D)
    @unpack gamma = equations
    rho = 2.0 + sin(pi * x[1] / 5.0) * sin(pi * x[2] / 5.0) * sin(pi * x[3] / 5.0)
    return prim2cons(SVector{5}(rho, .1, .2, .3, rho^gamma), equations)
end

run_benchmarks(run_benchmarks_2d, polydegs, flux_ranocha, initial_condition_sine, "sine")
run_benchmarks(run_benchmarks_3d, polydegs, flux_ranocha, initial_condition_sine, "sine")

## compute relative performance
using DelimitedFiles
using Measurements

function compute_relative_performance(filename)
  header = "# Polydeg Precomputed logarithms/primitive-variables-mean std"
  data = readdlm(filename, comments=true)
  cons = data[:, 2] .± data[:, 3]
  prim = data[:, 4] .± data[:, 5]
  prim_logs = data[:, 6] .± data[:, 7]  
  ratio = prim_logs ./ prim
  open(filename[1:end-4] * "_relative.dat", "w") do io
    println(io, header)
    writedlm(io, hcat(data[:, 1], # polydeg
                      Measurements.value.(ratio),
                      Measurements.uncertainty.(ratio)))
  end
end

compute_relative_performance(joinpath(@__DIR__, "precomputed_2D_flux_ranocha_vortex.dat"))
compute_relative_performance(joinpath(@__DIR__, "precomputed_3D_flux_ranocha_vortex.dat"))
compute_relative_performance(joinpath(@__DIR__, "precomputed_2D_flux_ranocha_random.dat"))
compute_relative_performance(joinpath(@__DIR__, "precomputed_3D_flux_ranocha_random.dat"))
compute_relative_performance(joinpath(@__DIR__, "precomputed_2D_flux_ranocha_sine.dat"))
compute_relative_performance(joinpath(@__DIR__, "precomputed_3D_flux_ranocha_sine.dat"))
