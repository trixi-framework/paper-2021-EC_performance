
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using DelimitedFiles, Statistics
using OrdinaryDiffEq
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


function run_benchmarks(polydeg, volume_flux, surface_flux, mesh)
  println("Benchmarking (", volume_flux, ", ", surface_flux, ") in ",
          ndims(mesh), "D with polydeg=", polydeg, " on ", nameof(typeof(mesh)))

  if ndims(mesh) == 2
    equations = CompressibleEulerEquations2D(1.4)
  else
    equations = CompressibleEulerEquations3D(1.4)
  end

  solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
                 volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  semi = SemidiscretizationHyperbolic(mesh, equations,
                                      initial_condition_isentropic_vortex, solver)

  # compile
  ode = semidiscretize(semi, (0.0, 1.0))
  solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
        save_everystep=false, maxiters=1)

  pids = zeros(5)
  for i in eachindex(pids)
    sleep(0.5)
    take!(semi.performance_counter)
    solve(ode, RDPK3SpFSAL49(), abstol=1.0e-8, reltol=1.0e-8,
          save_everystep=false, maxiters=50)
    pid = 1.0e-9 * take!(semi.performance_counter) / Trixi.ndofs(semi)
    pids[i] = pid
  end

  return mean(pids), std(pids)
end


function run_benchmarks(polydegs, ndims::Integer, volume_flux, surface_flux)
  header = "# Polydeg"
  header = header * "\tTreeMesh-mean\tstd"
  header = header * "\tStructuredMesh-mean\tstd"
  header = header * "\tP4estMesh-mean\tstd"
  filename = joinpath(@__DIR__, "pids_$(ndims)D_$(volume_flux)_$(first(polydegs))_$(last(polydegs)).dat")
  open(filename, "w") do io
    println(io, header)
  end

  if ndims == 2
    coordinates_min = (-5.0, -5.0)
    coordinates_max = ( 5.0,  5.0)
    initial_refinement_level = 3
    cells_per_dimension = 2^initial_refinement_level .* (1, 1)
  else
    coordinates_min = (-5.0, -5.0, -5.0)
    coordinates_max = ( 5.0,  5.0,  5.0)
    initial_refinement_level = 3
    cells_per_dimension = 2^initial_refinement_level .* (1, 1, 1)
  end

  for polydeg in polydegs
    results_line = Float64[polydeg]

    # TreeMesh
    let mesh = TreeMesh(coordinates_min, coordinates_max,
                        initial_refinement_level=initial_refinement_level,
                        n_cells_max=100_000, periodicity=true)

      pid_mean, pid_std = run_benchmarks(polydeg, volume_flux, surface_flux, mesh)
      push!(results_line, pid_mean)
      push!(results_line, pid_std)
    end

    # StructuredMesh
    let mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                              periodicity=true)

      pid_mean, pid_std = run_benchmarks(polydeg, volume_flux, surface_flux, mesh)
      push!(results_line, pid_mean)
      push!(results_line, pid_std)
    end

    # P4estMesh
    let mesh = P4estMesh(cells_per_dimension; coordinates_min, coordinates_max,
                         polydeg=1, periodicity=true)

      pid_mean, pid_std = run_benchmarks(polydeg, volume_flux, surface_flux, mesh)
      push!(results_line, pid_mean)
      push!(results_line, pid_std)
    end

    open(filename, "a") do io
      writedlm(io, results_line')
    end
  end

  return filename
end


## run benchmarks
polydegs = 3:15
# run_benchmarks(polydegs, 2, flux_shima_etal, flux_shima_etal)
# run_benchmarks(polydegs, 2, flux_ranocha, flux_ranocha)

# run_benchmarks(polydegs, 3, flux_shima_etal, flux_shima_etal)
# run_benchmarks(polydegs, 3, flux_ranocha, flux_ranocha)
