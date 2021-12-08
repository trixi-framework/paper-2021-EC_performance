
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using DelimitedFiles, LinearAlgebra, Statistics
using BenchmarkTools
using Trixi


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


# 2D version
function overintegration_kernel!(du_overinte, du_fluxdiff, u_overinte, u_fluxdiff,
                                 lo2hi, tmp_lo2hi, hi2lo, tmp_hi2lo,
                                 element, mesh,
                                 nonconservative_terms, equations,
                                 dg, cache)

  Trixi.multiply_dimensionwise!(view(u_overinte, :, :, :, element),
                                lo2hi,
                                view(u_fluxdiff, :, :, :, element),
                                tmp_lo2hi)

  Trixi.weak_form_kernel!(du_overinte, u_overinte, element, mesh,
                          nonconservative_terms, equations, dg, cache)

  Trixi.multiply_dimensionwise!(view(du_fluxdiff, :, :, :, element),
                                hi2lo,
                                view(du_overinte, :, :, :, element),
                                tmp_hi2lo)

  return nothing
end

function run_benchmarks_2d(polydeg_fluxdiff, polydeg_overinte, numerical_flux)
  println("Benchmarking ", numerical_flux, " in 2D with polydeg=", polydeg_fluxdiff)

  # setup basic code
  equations = CompressibleEulerEquations2D(1.4)

  coordinates_min = (-5.0, -5.0)
  coordinates_max = ( 5.0,  5.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=0, n_cells_max=100_000, periodicity=true)

  solver_fluxdiff = DGSEM(polydeg=polydeg_fluxdiff, surface_flux=numerical_flux,
                          volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))

  solver_overinte = DGSEM(polydeg=polydeg_overinte, surface_flux=numerical_flux,
                          volume_integral=VolumeIntegralWeakForm())

  semi_fluxdiff = SemidiscretizationHyperbolic(mesh, equations,
                                              initial_condition_isentropic_vortex, solver_fluxdiff)

  semi_overinte = SemidiscretizationHyperbolic(mesh, equations,
                                              initial_condition_isentropic_vortex, solver_overinte)

  u_ode_fluxdiff  = compute_coefficients(0.0, semi_fluxdiff)
  du_ode_fluxdiff = similar(u_ode_fluxdiff)
  u_ode_overinte  = compute_coefficients(0.0, semi_overinte)
  du_ode_overinte = similar(u_ode_overinte)
  GC.@preserve u_ode_fluxdiff du_ode_fluxdiff u_ode_overinte du_ode_overinte begin
    u_fluxdiff  = Trixi.wrap_array(u_ode_fluxdiff,  semi_fluxdiff)
    du_fluxdiff = Trixi.wrap_array(du_ode_fluxdiff, semi_fluxdiff)
    u_overinte  = Trixi.wrap_array(u_ode_overinte,  semi_overinte)
    du_overinte = Trixi.wrap_array(du_ode_overinte, semi_overinte)
    element = 1

    # benchmark flux differencing volume terms
    result_fluxdiff = @benchmark Trixi.split_form_kernel!(
      $du_fluxdiff, $u_fluxdiff, $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver_fluxdiff, $(semi_fluxdiff.cache)) setup=(fill!($du_fluxdiff, zero(eltype($du_fluxdiff))))
    display(result_fluxdiff)

    # benchmark overintegration volume terms
    nodes_lo = solver_fluxdiff.basis.nodes
    nodes_hi = solver_overinte.basis.nodes
    lo2hi = Trixi.polynomial_interpolation_matrix(nodes_lo, nodes_hi)

    V_lo, _ = Trixi.vandermonde_legendre(nodes_lo)
    V_lo = [V_lo zeros(length(nodes_lo), length(nodes_hi) - length(nodes_lo))]
    V_hi, _ = Trixi.vandermonde_legendre(nodes_hi)
    d = zeros(length(nodes_hi)); d[1:length(nodes_lo)] .= 1
    proj = Diagonal(d)
    hi2lo = V_lo * proj / V_hi

    tmp_lo2hi = zeros(size(u_fluxdiff, 1), size(lo2hi)...)
    tmp_hi2lo = zeros(size(u_fluxdiff, 1), size(hi2lo)...)

    result_overinte = @benchmark overintegration_kernel!(
      $du_overinte, $du_fluxdiff, $u_overinte, $u_fluxdiff,
      $lo2hi, $tmp_lo2hi, $hi2lo, $tmp_hi2lo,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $solver_overinte, $(semi_overinte.cache)) setup=(fill!($du_overinte, zero(eltype($du_overinte))))
    display(result_overinte)
  end
  println()

  return result_fluxdiff, result_overinte
end

# 3D version
function overintegration_kernel!(du_overinte, du_fluxdiff, u_overinte, u_fluxdiff,
                                 lo2hi, tmp1_lo2hi, tmp2_lo2hi,
                                 hi2lo, tmp1_hi2lo, tmp2_hi2lo,
                                 element, mesh,
                                 nonconservative_terms, equations,
                                 dg, cache)

  Trixi.multiply_dimensionwise!(view(u_overinte, :, :, :, :, element),
                                lo2hi,
                                view(u_fluxdiff, :, :, :, :, element),
                                tmp1_lo2hi, tmp2_lo2hi)

  Trixi.weak_form_kernel!(du_overinte, u_overinte, element, mesh,
                          nonconservative_terms, equations, dg, cache)

  Trixi.multiply_dimensionwise!(view(du_fluxdiff, :, :, :, :, element),
                                hi2lo,
                                view(du_overinte, :, :, :, :, element),
                                tmp1_hi2lo, tmp2_hi2lo)

  return nothing
end

function run_benchmarks_3d(polydeg_fluxdiff, polydeg_overinte, numerical_flux)
  println("Benchmarking ", numerical_flux, " in 3D with polydeg=", polydeg_fluxdiff)

  # setup basic code
  equations = CompressibleEulerEquations3D(1.4)

  coordinates_min = (-5.0, -5.0, -5.0)
  coordinates_max = ( 5.0,  5.0,  5.0)
  mesh = TreeMesh(coordinates_min, coordinates_max,
                  initial_refinement_level=0, n_cells_max=100_000, periodicity=true)

  solver_fluxdiff = DGSEM(polydeg=polydeg_fluxdiff, surface_flux=numerical_flux,
                          volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))

  solver_overinte = DGSEM(polydeg=polydeg_overinte, surface_flux=numerical_flux,
                          volume_integral=VolumeIntegralWeakForm())

  semi_fluxdiff = SemidiscretizationHyperbolic(mesh, equations,
                                              initial_condition_isentropic_vortex, solver_fluxdiff)

  semi_overinte = SemidiscretizationHyperbolic(mesh, equations,
                                              initial_condition_isentropic_vortex, solver_overinte)

  u_ode_fluxdiff  = compute_coefficients(0.0, semi_fluxdiff)
  du_ode_fluxdiff = similar(u_ode_fluxdiff)
  u_ode_overinte  = compute_coefficients(0.0, semi_overinte)
  du_ode_overinte = similar(u_ode_overinte)
  GC.@preserve u_ode_fluxdiff du_ode_fluxdiff u_ode_overinte du_ode_overinte begin
    u_fluxdiff  = Trixi.wrap_array(u_ode_fluxdiff,  semi_fluxdiff)
    du_fluxdiff = Trixi.wrap_array(du_ode_fluxdiff, semi_fluxdiff)
    u_overinte  = Trixi.wrap_array(u_ode_overinte,  semi_overinte)
    du_overinte = Trixi.wrap_array(du_ode_overinte, semi_overinte)
    element = 1

    # benchmark flux differencing volume terms
    result_fluxdiff = @benchmark Trixi.split_form_kernel!(
      $du_fluxdiff, $u_fluxdiff, $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $numerical_flux, $solver_fluxdiff, $(semi_fluxdiff.cache)) setup=(fill!($du_fluxdiff, zero(eltype($du_fluxdiff))))
    display(result_fluxdiff)

    # benchmark overintegration volume terms
    nodes_lo = solver_fluxdiff.basis.nodes
    nodes_hi = solver_overinte.basis.nodes
    lo2hi = Trixi.polynomial_interpolation_matrix(nodes_lo, nodes_hi)

    V_lo, _ = Trixi.vandermonde_legendre(nodes_lo)
    V_lo = [V_lo zeros(length(nodes_lo), length(nodes_hi) - length(nodes_lo))]
    V_hi, _ = Trixi.vandermonde_legendre(nodes_hi)
    d = zeros(length(nodes_hi)); d[1:length(nodes_lo)] .= 1
    proj = Diagonal(d)
    hi2lo = V_lo * proj / V_hi

    tmp1_lo2hi = zeros(size(u_fluxdiff, 1), size(lo2hi, 1), size(lo2hi, 2), size(lo2hi, 2))
    tmp2_lo2hi = zeros(size(u_fluxdiff, 1), size(lo2hi, 1), size(lo2hi, 1), size(lo2hi, 2))
    tmp1_hi2lo = zeros(size(u_fluxdiff, 1), size(hi2lo, 1), size(hi2lo, 2), size(hi2lo, 2))
    tmp2_hi2lo = zeros(size(u_fluxdiff, 1), size(hi2lo, 1), size(hi2lo, 1), size(hi2lo, 2))

    result_overinte = @benchmark overintegration_kernel!(
      $du_overinte, $du_fluxdiff, $u_overinte, $u_fluxdiff,
      $lo2hi, $tmp1_lo2hi, $tmp2_lo2hi, $hi2lo, $tmp1_hi2lo, $tmp2_hi2lo,
      $element, $mesh,
      $(Trixi.have_nonconservative_terms(equations)), $equations,
      $solver_overinte, $(semi_overinte.cache)) setup=(fill!($du_overinte, zero(eltype($du_overinte))))
    display(result_overinte)
  end
  println()

  return result_fluxdiff, result_overinte
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
  means = zeros(length(polydegs))
  stds  = zero(means)
  for (i, polydeg) in enumerate(polydegs)
    result, _ = run_benchmarks_inner(polydeg, polydeg, numerical_flux)
    means[i] = time(mean(result)) * 1.0e-9
    stds[i]  = time(std(result))  * 1.0e-9
  end
  results = hcat(results, means, stds)
  header = header * "\tFluxDifferencing-mean\tstd"

  # overintegration, polydeg + 1
  println("#"^80)
  println("overintegration, polydeg + 1\n")
  means = zeros(length(polydegs))
  stds  = zero(means)
  for (i, polydeg) in enumerate(polydegs)
    _, result = run_benchmarks_inner(polydeg, polydeg + 1, flux_central)
    means[i] = time(mean(result)) * 1.0e-9
    stds[i]  = time(std(result))  * 1.0e-9
  end
  results = hcat(results, means, stds)
  header = header * "\tpolydeg+1-mean\tstd"

  # overintegration, polydeg + 2
  println("#"^80)
  println("overintegration, polydeg + 2\n")
  means = zeros(length(polydegs))
  stds  = zero(means)
  for (i, polydeg) in enumerate(polydegs)
    _, result = run_benchmarks_inner(polydeg, polydeg + 2, flux_central)
    means[i] = time(mean(result)) * 1.0e-9
    stds[i]  = time(std(result))  * 1.0e-9
  end
  results = hcat(results, means, stds)
  header = header * "\tpolydeg+2-mean\tstd"

  # overintegration, 3 * polydeg ÷ 2
  println("#"^80)
  println("overintegration, 3 * polydeg ÷ 2\n")
  means = zeros(length(polydegs))
  stds  = zero(means)
  for (i, polydeg) in enumerate(polydegs)
    _, result = run_benchmarks_inner(polydeg, (3 * polydeg) ÷ 2, flux_central)
    means[i] = time(mean(result)) * 1.0e-9
    stds[i]  = time(std(result))  * 1.0e-9
  end
  results = hcat(results, means, stds)
  header = header * "\t3*polydeg÷2-mean\tstd"

  # overintegration, 2 * polydeg
  println("#"^80)
  println("overintegration, 2 * polydeg\n")
  means = zeros(length(polydegs))
  stds  = zero(means)
  for (i, polydeg) in enumerate(polydegs)
    _, result = run_benchmarks_inner(polydeg, 2 * polydeg, flux_central)
    means[i] = time(mean(result)) * 1.0e-9
    stds[i]  = time(std(result))  * 1.0e-9
  end
  results = hcat(results, means, stds)
  header = header * "\t2*polydeg-mean\tstd"

  # write results to a file
  open(joinpath(@__DIR__, "overintegration_$(ndims(run_benchmarks_inner))D_$(numerical_flux).dat"), "w") do io
    println(io, header)
    writedlm(io, results)
  end

  return results
end

function run_benchmarks_2d(polydegs, numerical_flux)
  run_benchmarks(run_benchmarks_2d, polydegs, numerical_flux)
end

function run_benchmarks_3d(polydegs, numerical_flux)
  run_benchmarks(run_benchmarks_3d, polydegs, numerical_flux)
end



## run benchmarks
polydegs = 3:15
run_benchmarks_2d(polydegs, flux_shima_etal)
run_benchmarks_2d(polydegs, flux_ranocha)

run_benchmarks_3d(polydegs, flux_shima_etal)
run_benchmarks_3d(polydegs, flux_ranocha)
