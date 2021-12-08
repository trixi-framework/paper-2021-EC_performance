# Running EC performance tests with FLUXO

## Requirements
* Fortran compiler
* HDF5 with Fortran support
* CMake (minimum required version is 3.5.1)
* Python (>=3) with numpy

## Instructions

To run the performance tests, follow the instructions:

* Move to this directory and clone the FLUXO repository:
  ```
  cd 2021_EC_performance/code/ec_performance_fluxo
  git clone git@github.com:project-fluxo/fluxo.git
  ```
* Move to the branch of FLUXO were the performance improvements are implemented:
  ```
  cd fluxo
  git checkout performance
  cd ..
  ```
* *Optional:* By default, FLUXO will compile with the Intel compiler using the
  `-xHost` flag, i.e., it will try to generate optimal machine code for the CPU
  used during compile time. If you want to specifically compile for a given
  architecture or restrict to a specific instruction set, you can change this in
  FLUXO's `CMakeLists.txt` file
  by modifying the
  [relevant lines](https://github.com/project-fluxo/fluxo/blob/01b8bfc62a566dde12769989c87630754deab7e1/CMakeLists.txt#L389-L439)
  that set the optimization flags for the respectively used compiler.
  For example, to make `ifort` only use `AVX2` instruction on an
  `AVX512`-enabled CPU, change the line
  ```cmake
        SET(FLUXO_INSTRUCTION "-xHost")
  ```
  to
  ```cmake
        SET(FLUXO_INSTRUCTION "-xCore-AVX2")
  ```
* *Optional:* For the tests, you can use the Cartesian mesh with 8² elements under `MeshGeneration/CartBoxPeriodic_008_008_008_mesh.h5`. If you prefer to generate the mesh, you can do it using [hopr](https://gitlab.com/project-fluxo/hopr):
  ```
  cd MeshGeneration
  hopr parameter_hopr.ini
  cd ..
  ```
* Go to the folder `serialTests` and run the tests using the python script:
  ```
  cd serialTests
  python3 runTests.py [options]
  ```
  Following optional arguments are available:
  ```
  -h, --help           show this help message and exit
  --compiler COMPILER  Fortran compiler for FLUXO. If empty, the Fortran compiler found by cmake is used.
  --id ID              Test Identifier. Used for the build directories ("build/build_Euler_<id>...") and the results files ("serialTests/results/results_<id>_.dat")
  --flux FLUX          Volume and surface flux used for the tests.
                           32: [DEFAULT] Ranocha flux
                           40: Shima et al. flux
  --precompileN        [DEFAULT] Compile fluxo with fixed polynomial degree.
  --no-precompileN     [optional] deactivate --precompileN.
  --precompileFlux     [DEFAULT] Compile fluxo with volume flux.
  --no-precompileFlux  [optional] deactivate --precompileN.
  --extra_options EXTRA_OPTIONS
                        Additional options for cmake. For example, to compile with Gauss nodes use --extra_options=" -DFLUXO_DISC_NODETYPE=GAUSS"

  ```
  
  In particular, the tests for the EC performance can be launched with the following commands:
  ```
  python3 runTests.py  --compiler ifort --id Ranocha --flux 32
  python3 runTests.py  --compiler ifort --id Ranocha_noN_noFlux --flux 32 --no-precompileN --no-precompileFlux
  python3 runTests.py  --compiler ifort --id Ranocha_noN --flux 32 --no-precompileN
  python3 runTests.py  --compiler ifort --id Ranocha_noFlux --flux 32 --no-precompileFlux
  python3 runTests.py  --compiler ifort --id Shima --flux 40
  python3 runTests.py  --compiler ifort --id Ranocha_Gauss --flux 32 --extra_options=" -DFLUXO_DISC_NODETYPE=GAUSS"
  ```
    
All results were obtained with [this version of fluxo](https://github.com/project-fluxo/fluxo/tree/01b8bfc62a566dde12769989c87630754deab7e1).
