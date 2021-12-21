# Efficient implementation of modern entropy stable and kinetic energy preserving discontinuous Galerkin methods for conservation laws

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5792576.svg)](https://doi.org/10.5281/zenodo.5792576)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{ranocha2021efficient,
  title={Efficient implementation of modern entropy stable and kinetic energy
         preserving discontinuous {G}alerkin methods for conservation laws},
  author={Ranocha, Hendrik and Schlottke-Lakemper, Michael and Chan, Jesse and
          Rueda-Ram\'{i}rez, Andr{\'e}s M and Winters, Andrew R and
          Hindenlang, Florian and Gassner, Gregor J},
  year={2021},
  month={12},
  eprint={2112.10517},
  eprinttype={arxiv},
  eprintclass={cs.MS}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{ranocha2021efficientRepro,
  title={Reproducibility repository for
         {E}fficient implementation of modern entropy stable and kinetic energy
         preserving discontinuous {G}alerkin methods for conservation laws},
  author={Ranocha, Hendrik and Schlottke-Lakemper, Michael and Chan, Jesse and
          Rueda-Ram\'{i}rez, Andr{\'e}s M and Winters, Andrew R and
          Hindenlang, Florian and Gassner, Gregor J},
  year={2021},
  month={12},
  howpublished={\url{https://github.com/trixi-framework/paper-2021-EC\_performance}},
  doi={10.5281/zenodo.5792576}
}
```


## Abstract

Many modern discontinuous Galerkin (DG) methods for conservation laws make use of summation by parts operators and flux differencing to achieve kinetic energy preservation or entropy stability. While these techniques increase the robustness of DG methods significantly, they are also computationally more demanding than standard weak form nodal DG methods. We present several implementation techniques to improve the efficiency of flux differencing DG methods that use tensor product quadrilateral or hexahedral elements, in 2D or 3D respectively. Focus is mostly given to CPUs and DG methods for the compressible Euler equations, although these techniques are generally also useful for GPU computing and other physical systems including the compressible Navier-Stokes and magnetohydrodynamics equations. We present results using two open source codes, Trixi.jl written in Julia and FLUXO written in Fortran, to demonstrate that our proposed implementation techniques are applicable to different code bases and programming languages.


## Numerical experiments

The numerical experiments presented in the paper use [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
and [FLUXO](https://gitlab.com/project-fluxo/fluxo).
To reproduce the numerical experiments using Trixi.jl, you need to install
[Julia](https://julialang.org/).

The subfolders of this repository contain `README.md` files with instructions
to reproduce the numerical experiments, including postprocessing.

The numerical experiments were carried out using Julia v1.7.0.


## Authors

- [Hendrik Ranocha](https://ranocha.de) (University of Münster, Germany)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (University of Stuttgart, Germany)
- [Jesse Chan](https://jlchan.github.io) (Rice University, USA)
- [Andrés M. Rueda-Ramírez](https://www.mi.uni-koeln.de/NumSim/dr-andres-rueda-ramirez) (University of Cologne, Germany)
- [Andrew Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden)
- Florian Hindenlang (Max Planck Institute for Plasma Physics, Germany)
- [Gregor Gassner](https://www.mi.uni-koeln.de/NumSim/gregor-gassner) (University of Cologne, Germany)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
