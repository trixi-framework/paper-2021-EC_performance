
# Gather results from the Vulcan cluster at HLRS.
# The setup is described in `vulcan.md`. This script postprocesses
# the data generated that way by moving output files to get the
# same structure as obtained by running everything locally as
# described in the `README.md` files in each subdirectory.

using DelimitedFiles

function postprocess_vulcan(directory)
  for (root, dirs, files) in walkdir(directory)
    for dir in dirs
      # directories are named "run_IDENTIFIER"
      if !startswith(dir, "run_")
        continue
      end

      id = dir[5:end]
      jobscript = "jobscript_" * id * ".sh"
      if !(jobscript in files) && !("jobscript.sh" in files)
        @warn "There is no jobscript '$jobscript' associated to '$dir' in $root"
        continue
      end

      for (subroot, subdirs, subfiles) in walkdir(joinpath(root, dir))
        for file in subfiles
          if endswith(file, ".dat") || endswith(file, ".txt")
            src = joinpath(subroot, file)
            dst = joinpath(root, file)
            @info "Copying $src to $dst"
            cp(src, dst, force=true)
          end
        end
      end
    end
  end
end


postprocess_vulcan(joinpath(@__DIR__, "numerical_fluxes"))
postprocess_vulcan(joinpath(@__DIR__, "overintegration"))
postprocess_vulcan(joinpath(@__DIR__, "precompute_logs"))
postprocess_vulcan(joinpath(@__DIR__, "primitive_variables"))
