
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


function run_benchmarks(equations, polydeg, volume_flux, surface_flux, approximation_type)
  println("Benchmarking (", volume_flux, ", ", surface_flux, ") in ",
          ndims(equations), "D with polydeg=", polydeg, " and ", approximation_type)

  if ndims(equations) == 2
    element_type = Quad()
  else
    element_type = Hex()
  end
  solver = DGMulti(polydeg=polydeg, element_type=element_type,
                   approximation_type=approximation_type,
                   surface_integral=SurfaceIntegralWeakForm(surface_flux),
                   volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

  vertex_coordinates, EToV = StartUpDG.uniform_mesh(element_type, 8)
  vertex_coordinates = map(x -> 5 .* x, vertex_coordinates) # map domain to [-5, 5]^ndims(equations)
  mesh = VertexMappedMesh(vertex_coordinates, EToV, solver,
                          is_periodic=ntuple(_ -> true, ndims(equations)))

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


function run_benchmarks(polydegs, ndims, volume_flux, surface_flux)
  header = "# Polydeg"
  header = header * "\tSBP-mean\tstd"
  header = header * "\tGaussSBP-mean\tstd"
  filename = joinpath(@__DIR__, "Gauss_$(ndims)D_$(volume_flux)_$(first(polydegs))_$(last(polydegs)).dat")
  open(filename, "w") do io
    println(io, header)
  end

  if ndims == 2
    equations = CompressibleEulerEquations2D(1.4)
  else
    equations = CompressibleEulerEquations3D(1.4)
  end

  for polydeg in polydegs
    results_line = Float64[polydeg]

    # SBP
    let approximation_type = SBP()
      pid_mean, pid_std = run_benchmarks(equations, polydeg, volume_flux, surface_flux,
                                         approximation_type)
      push!(results_line, pid_mean)
      push!(results_line, pid_std)
    end

    # GaussSBP
    let approximation_type = GaussSBP()
      pid_mean, pid_std = run_benchmarks(equations, polydeg, volume_flux, surface_flux,
                                         approximation_type)
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
# run_benchmarks(polydegs, 2, flux_ranocha, flux_ranocha)
# run_benchmarks(polydegs, 3, flux_ranocha, flux_ranocha)
