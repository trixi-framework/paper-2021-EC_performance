
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using TOML
using Printf

using LIKWID
using Trixi

using PyCall
import PyPlot; plt = PyPlot

cycler = pyimport("cycler").cycler
line_cycler   = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]))
marker_cycler = (cycler(color=["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]) +
                 cycler(linestyle=["none", "none", "none", "none", "none", "none", "none"]) +
                 cycler(marker=["4", "2", "3", "1", "+", "x", "."]))

plt.rc("axes", prop_cycle=line_cycler)
# plt.rc("text", usetex=true)
# plt.rc("text.latex", preamble="\\usepackage{newpxtext}\\usepackage{newpxmath}\\usepackage{commath}\\usepackage{mathtools}")
plt.rc("font", family="serif", size=18.)
plt.rc("savefig", dpi=100)
plt.rc("legend", loc="best", fontsize="medium", fancybox=true, framealpha=0.5)
plt.rc("lines", linewidth=2.5, markersize=10, markeredgewidth=2.5)


## setup data dictionary to save results
data = Dict{String, Any}()


## gather data for the empirical roofline model

# measure optimistic peakflops (AVX2 FMA or AVX512 FMA if available)
L1_cache_size = LIKWID.get_cpu_topology().cacheLevels[1].size ÷ 1024 # in kB

cpuinfo = LIKWID.get_cpu_info()
if occursin("AVX512", cpuinfo.features)
  likwid_bench_kernel_flops = "peakflops_avx512_fma"
elseif occursin("AVX2", cpuinfo.features)
  likwid_bench_kernel_flops = "peakflops_avx_fma"
else
  likwid_bench_kernel_flops = "peakflops_sse_fma"
end

max_flops_string = read(`likwid-bench -t $likwid_bench_kernel_flops -W N:$(L1_cache_size)kB:1`, String)
max_flops = parse(Float64, match(r"(MFlops/s:\s+)(\d+\.\d+)", max_flops_string).captures[2]) / 1024
data["max_flops"] = max_flops

# measure optimistic memory bandwidth using reads
if occursin("AVX512", cpuinfo.features)
  likwid_bench_kernel_bandwidth = "load_avx512"
elseif occursin("AVX2", cpuinfo.features)
  likwid_bench_kernel_bandwidth = "load_avx"
else
  likwid_bench_kernel_bandwidth = "load_sse"
end

max_bandwidth_string = read(`likwid-bench -t $likwid_bench_kernel_bandwidth -W N:2GB:1`, String)
max_bandwidth = parse(Float64, match(r"(MByte/s:\s+)(\d+\.\d+)", max_bandwidth_string).captures[2]) / 1024
data["max_bandwidth"] = max_bandwidth


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

function intensity_flops(numerical_flux, ::Type{MeshType}) where {MeshType}
  equations = CompressibleEulerEquations3D(1.4)
  solver = DGSEM(polydeg=3, surface_flux=numerical_flux,
                 volume_integral=VolumeIntegralFluxDifferencing(numerical_flux))
  coordinates_min = (-5.0, -5.0, -5.0)
  coordinates_max = ( 5.0,  5.0,  5.0)
  initial_refinement_level = 3
  cells_per_dimension = 2^initial_refinement_level .* (1, 1, 1)
  t = 0.0

  if MeshType == TreeMesh
    mesh = TreeMesh(coordinates_min, coordinates_max,
                    initial_refinement_level=initial_refinement_level,
                    n_cells_max=100_000, periodicity=true)
  elseif MeshType == StructuredMesh
    mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
                          periodicity=true)
  elseif MeshType == P4estMesh
    mesh = P4estMesh(cells_per_dimension; coordinates_min, coordinates_max,
                     polydeg=1, periodicity=true)
  else
    throw(ArgumentError("Mesh type $MeshType not supported."))
  end

  semi = SemidiscretizationHyperbolic(mesh, equations,
                                      initial_condition_isentropic_vortex, solver)

  u_ode = Trixi.compute_coefficients(t, semi)
  du_ode = zero(u_ode)

  GC.@preserve u_ode du_ode begin
    u = Trixi.wrap_array(u_ode, semi)
    du = Trixi.wrap_array(du_ode, semi)

    # warmup
    many_volume_terms!(du, u, semi, t, 1)
    sleep(0.5)

    metrics, events = @perfmon "MEM_DP" many_volume_terms!(du, u, semi, t, 1)
    intensity = metrics["MEM_DP"][1]["Operational intensity"]
    sleep(0.5)

    metrics, events = @perfmon "MEM_DP" many_volume_terms!(du, u, semi, t, 5000)
    flops = metrics["MEM_DP"][1]["DP [MFLOP/s]"] / 1024
  end

  return intensity, flops
end


## gather data for volume terms implemented in Trixi.jl
fluxes = (flux_shima_etal, flux_shima_etal_turbo, flux_ranocha, flux_ranocha_turbo)
meshes = (TreeMesh, StructuredMesh, P4estMesh)

counter = 1
for numerical_flux in fluxes
  for mesh in meshes
    intensity, flops = intensity_flops(numerical_flux, mesh)
    @info "Results" counter numerical_flux mesh intensity flops
    new_data = Dict{String,Any}()
    new_data["numerical_flux"] = string(numerical_flux)
    new_data["mesh"] = string(mesh)
    new_data["intensity"] = intensity
    new_data["flops"] = flops
    data[@sprintf("Experiment %02d", counter)] = new_data
    global counter += 1
  end
end

open("results.txt", "w") do io
  TOML.print(io, data; sorted=true)
end

