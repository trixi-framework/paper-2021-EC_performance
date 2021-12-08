# Running benchmarks

Documentation for how benchmarks were run on HLRS' Vulcan cluster.

## Job submission
All jobscripts are named `jobscript_...sh` and should be submitted using `qsub`
from the folder in which they are located. For example, for `jobscript.sh`,
execute
```bash
qsub jobscript.sh
```
The Trixi.jl tests can all be submitted simultaneously, while the FLUXO results
need to be submitted serially (or a separate `serialTests` directory needs to be
created for each run).

## Trixi.jl

To pre-instantiate all benchmark directories (since we have no internet
connection on the compute nodes), run the following command first:
```bash
 julia -e 'include("preinstantiate.jl"); preinstantiate()'
```
All resulting files for each run were gathered in the corresponding `run_XXX`
folder. Previous result files have **not been modified**.

For all tests, the pre-built binaries for Julia 1.7.0 were used.

To postprocess results from Vulcan and obtain the same directory structure as
from a local run following the instructions in the individual README.md files,
execute
```bash
julia postprocess_results_from_vulcan.jl
```

*Note:* For `Gauss` and `inlining`, the 3D/`flux_ranocha` runs for polydegs
3-13 did not complete (in `Gauss` there was a bus error, in `inlining` it exceeded the
requested walltime). In both cases, the results for polydeg 12 were still
written to the results file. Thus, this was renamed `part1a` in both cases, and
a separate run for polydeg 13 was started, named `part1b`.

**Warning:** For `Gauss`|3D|`flux_ranocha`|polydeg=13 *no* results have been obtained yet,
since the simulation consistently crashes after ~38 minutes with a `Bus error`.
I have not yet been able to figure out why :-/. Full error message:
```
signal (7): Bus error
in expression starting at none:1
macro expansion at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/VectorizationBase/OEl8L/src/llvm_intrin/vbroadcast.jl:74 [inlined]
_vbroadcast at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/VectorizationBase/OEl8L/src/llvm_intrin/vbroadcast.jl:80 [inlined]
vbroadcast at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/VectorizationBase/OEl8L/src/llvm_intrin/vbroadcast.jl:96 [inlined]
macro expansion at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/LoopVectorization/kVenK/src/reconstruct_loopset.jl:713 [inlined]
_turbo_! at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/LoopVectorization/kVenK/src/reconstruct_loopset.jl:713 [inlined]
ploopmul! at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/Octavian/Rlmrt/src/macrokernels.jl:30 [inlined]
packaloopmul! at /zhome/academic/HLRS/hlrs/hpcschlo/.julia/packages/Octavian/Rlmrt/src/macrokernels.jl:138
unknown function (ip: 0xbeaacd687bdba98b)
Allocations: 168505022 (Pool: 168391374; Big: 113648); GC: 84
/var/spool/PBS/mom_priv/jobs/466561.cl5intern.SC: line 9: 61945 Bus error               julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(13, 3, flux_ranocha, flux_ranocha)'

real    38m50.716s
user    38m32.338s
sys     0m4.634s
```

## FLUXO
All resulting *auxiliary* files for each run were gathered in the corresponding
`run_XXX` folder in the `serialTests` directory. Previous  result files **have
been overwritten**.

The following configuration was used for FLUXO:
* Source: https://github.com/project-fluxo/fluxo
* Branch: `performance`
* Commit hash: `01b8bfc62a566dde12769989c87630754deab7e1`
* Intel Fortran Compiler 19.1.3
