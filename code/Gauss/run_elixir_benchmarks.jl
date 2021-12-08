
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using Trixi, OrdinaryDiffEq
using Trixi: TimerOutputs
using DelimitedFiles


## set up benchmark code
function compute_relative_performance(to, semi)
    pid(to, semi) = 1e-9 * TimerOutputs.time(to) / (TimerOutputs.ncalls(to) * Trixi.ndofs(semi))

    time_total = pid(to["rhs!"], semi)
    time_volume = pid(to["rhs!"]["volume integral"], semi)
    time_entropy_projection = pid(to["rhs!"]["entropy_projection!"], semi)
    time_surface = pid(to["rhs!"]["surface integral"], semi)
    time_interface = pid(to["rhs!"]["interface flux"], semi)
    time_other = time_total - (time_volume + time_entropy_projection +
                               time_interface + time_surface)

    open(joinpath(@__DIR__, "Gauss_$(ndims(semi))D_relative_timings.dat"), "a") do io
        writedlm(io, hcat(Trixi.polydeg(semi.solver), # polydeg
                          time_volume,
                          time_entropy_projection,
                          time_surface,
                          time_interface,
                          time_other)
                )
    end
end


function generate_timings(elixir_path, polydeg_range)
    # run for a very short time to avoid timing precompilation
    redirect_stdout(devnull) do
        trixi_include(elixir_path, polydeg=first(polydeg_range), tspan=(0.0, 1e-5))
    end

    for polydeg in polydeg_range
        # run for a very short time to avoid timing precompilation
        redirect_stdout(devnull) do
            trixi_include(elixir_path, polydeg=polydeg, tspan=(0.0, 1e-5))
        end

        t = @elapsed begin
            redirect_stdout(devnull) do
                trixi_include(elixir_path, polydeg=polydeg)
            end
        end
        println("polydeg=$polydeg completed in $t seconds")

        if polydeg==first(polydeg_range)
            open(joinpath(@__DIR__, "Gauss_$(ndims(semi))D_relative_timings.dat"), "w") do io
                header = "Polydeg\ttime_volume\ttime_entropy_projection\ttime_surface\ttime_interface\tother"
                println(io, header)
            end
        end
        to = Trixi.timer()
        compute_relative_performance(to, semi)
    end
end

generate_timings(joinpath(@__DIR__, "elixir_khi_2D.jl"), 3:15)
generate_timings(joinpath(@__DIR__, "elixir_tgv_3D.jl"), 3:15)
