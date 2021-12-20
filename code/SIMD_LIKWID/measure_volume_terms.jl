## Run this file using LIKWID:
# ```bash
# likwid-perfctr -C 0 -g MEM_DP -m julia --check-bounds=no --threads=1 measure_volume_terms.jl
# ```
# You can also inspect other groups, e.g., `-g FLOPS_DP`.


## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using LIKWID
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


function many_volume_terms!(du, u, semi, t, n)
  nonconservative_terms = Trixi.have_nonconservative_terms(semi.equations)
  for _ in 1:n
    Trixi.calc_volume_integral!(du, u, semi.mesh,
      nonconservative_terms, semi.equations,
      semi.solver.volume_integral, semi.solver, semi.cache)
  end
end

function run_measurements(polydeg=3, n=5*10^3)

  equations = CompressibleEulerEquations3D(1.4)

  Marker.init()

  for numerical_flux in (flux_shima_etal_turbo, flux_ranocha_turbo)

    solver = DGSEM(polydeg=polydeg, surface_flux=numerical_flux,
                   volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))

    coordinates_min = (-5.0, -5.0, -5.0)
    coordinates_max = ( 5.0,  5.0,  5.0)
    initial_refinement_level = 3
    cells_per_dimension = 2^initial_refinement_level .* (1, 1, 1)
    t = 0.0


    # TreeMesh
    let
      mesh = TreeMesh(coordinates_min, coordinates_max,
                      initial_refinement_level=initial_refinement_level,
                      n_cells_max=100_000, periodicity=true)

      semi = SemidiscretizationHyperbolic(mesh, equations,
                                          initial_condition_isentropic_vortex, solver)

      u_ode = Trixi.compute_coefficients(t, semi)
      du_ode = zero(u_ode)

      GC.@preserve u_ode du_ode begin
        u = Trixi.wrap_array(u_ode, semi)
        du = Trixi.wrap_array(du_ode, semi)

        # compile and cool down
        many_volume_terms!(du, u, semi, t, 1)
        sleep(1.0)

        @region "TreeMesh-$(numerical_flux)" begin
          many_volume_terms!(du, u, semi, t, n)
        end
      end
    end


    # StructuredMesh
    let
      mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                            periodicity=true)

      semi = SemidiscretizationHyperbolic(mesh, equations,
                                          initial_condition_isentropic_vortex, solver)

      u_ode = Trixi.compute_coefficients(t, semi)
      du_ode = zero(u_ode)

      GC.@preserve u_ode du_ode begin
        u = Trixi.wrap_array(u_ode, semi)
        du = Trixi.wrap_array(du_ode, semi)

        # compile and cool down
        many_volume_terms!(du, u, semi, t, 1)
        sleep(1.0)

        @region "StructuredMesh-$(numerical_flux)" begin
          many_volume_terms!(du, u, semi, t, n)
        end
      end
    end


    # P4estMesh
    let
      mesh = P4estMesh(cells_per_dimension; coordinates_min, coordinates_max,
                       polydeg=1, periodicity=true)

      semi = SemidiscretizationHyperbolic(mesh, equations,
                                          initial_condition_isentropic_vortex, solver)

      u_ode = Trixi.compute_coefficients(t, semi)
      du_ode = zero(u_ode)

      GC.@preserve u_ode du_ode begin
        u = Trixi.wrap_array(u_ode, semi)
        du = Trixi.wrap_array(du_ode, semi)

        # compile and cool down
        many_volume_terms!(du, u, semi, t, 1)
        sleep(1.0)

        @region "P4estMesh-$(numerical_flux)" begin
          many_volume_terms!(du, u, semi, t, n)
        end
      end
    end

  end

  Marker.close()

  return nothing
end


run_measurements()
