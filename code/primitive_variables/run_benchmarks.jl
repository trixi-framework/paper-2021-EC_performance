
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using DelimitedFiles, LinearAlgebra, Statistics
using BenchmarkTools
using StaticArrays
using StrideArrays
using Trixi


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


function split_form_prim_kernel!(du, u_prim, u_cons,
                                 element, mesh,
                                 nonconservative_terms, equations,
                                 volume_flux::VF, dg, cache) where {VF}

  cons2prim!(vec(u_prim), reshape(u_cons, nvariables(equations), :),
             equations, dg)

  Trixi.split_form_kernel!(du, u_prim, element, mesh,
                           nonconservative_terms, equations,
                           volume_flux, dg, cache)

  return nothing
end


# 2D version
function run_benchmarks_2d(polydeg, numerical_flux)
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
                                      initial_condition_isentropic_vortex, solver)

  u_ode  = compute_coefficients(0.0, semi)
  du_ode = similar(u_ode)
  GC.@preserve u_ode du_ode begin
    u  = Trixi.wrap_array(u_ode,  semi)
    du = Trixi.wrap_array(du_ode, semi)
    u_prim = StrideArray(undef,
                         PrimitiveVariablesCompressibleEuler2D{eltype(u)},
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
  end
  println()

  return result_cons, result_prim
end


# 3D version
function run_benchmarks_3d(polydeg, numerical_flux)
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
                                      initial_condition_isentropic_vortex, solver)

  u_ode  = compute_coefficients(0.0, semi)
  du_ode = similar(u_ode)
  GC.@preserve u_ode du_ode begin
    u  = Trixi.wrap_array(u_ode,  semi)
    du = Trixi.wrap_array(du_ode, semi)
    u_prim = StrideArray(undef,
                         PrimitiveVariablesCompressibleEuler3D{eltype(u)},
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
  end
  println()

  return result_cons, result_prim
end


# general code
Base.ndims(::typeof(run_benchmarks_2d)) = 2
Base.ndims(::typeof(run_benchmarks_3d)) = 3

function run_benchmarks(run_benchmarks_inner, polydegs, numerical_flux)
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
  for (i, polydeg) in enumerate(polydegs)
    result_cons, result_prim = run_benchmarks_inner(polydeg, numerical_flux)
    means_cons[i] = time(mean(result_cons)) * 1.0e-9
    stds_cons[i]  = time(std(result_cons))  * 1.0e-9
    means_prim[i] = time(mean(result_prim)) * 1.0e-9
    stds_prim[i]  = time(std(result_prim))  * 1.0e-9
  end
  results = hcat(results, means_cons, stds_cons, means_prim, stds_prim)
  header = header * "\tConservative-variables-mean\tstd" *
                    "\tPrimitive-variables-mean\tstd"

  # write results to a file
  open(joinpath(@__DIR__, "primitive_$(ndims(run_benchmarks_inner))D_$(numerical_flux).dat"), "w") do io
    println(io, header)
    writedlm(io, results)
  end

  return results
end



## run benchmarks
polydegs = 3:15
run_benchmarks(run_benchmarks_2d, polydegs, flux_shima_etal)
run_benchmarks(run_benchmarks_2d, polydegs, flux_ranocha)

run_benchmarks(run_benchmarks_3d, polydegs, flux_shima_etal)
run_benchmarks(run_benchmarks_3d, polydegs, flux_ranocha)
