[![DOI](https://zenodo.org/badge/1114221142.svg)](https://doi.org/10.5281/zenodo.19199177)
[![arXiv](https://img.shields.io/badge/arXiv-2603.26212-b31b1b.svg)](https://arxiv.org/abs/2603.26212)


# PaperUnfittedFEMDarcy

This repository holds the Julia software and code to reproduce the results reported in Section 8.1 of the paper _"Divergence-free unfitted finite element discretisations for the Darcy problem"_ by _Santiago Badia, Anne Boschman, Alberto F. Martín, Erik Nilsson, Ricardo Ruiz-Baier, and Sara Zahedi_. In particular, it contains:

* The h-convergence test on a cut square in Section 8.1.1. [here](https://github.com/amboschman/PaperUnfittedFEMDarcy.jl/blob/main/test/FIG8-1_h-convergence-on-a-cut-square.jl);
* The sensitivity analysis with respect to the penalty parameter in Section 8.1.2. [here](https://github.com/amboschman/PaperUnfittedFEMDarcy.jl/blob/main/test/FIG8-2_sensitivity-with-respect-to-gamma.jl);
* The pressure robustness test in Section 8.1.3. [here](https://github.com/amboschman/PaperUnfittedFEMDarcy.jl/blob/main/test/FIG8-3_pressure-robustness.jl);
* The conditioning analysis with respect to cut-cell length in Section 8.1.4. [here](https://github.com/amboschman/PaperUnfittedFEMDarcy.jl/blob/main/test/FIG8-4_conditioning-with-respect-to-cut-cell-length.jl).

**IMPORTANT NOTE** 
To generate the full results presented in the paper, the parameter settings need to be changed as outlined in the comments of each test file. A selective range of parameters is chosen as the default so that the code can run efficiently on a system with more limited computational resources. 

## Installation instructions

To run the code, the user has to generate a Manifest.toml file in such way that the code can be run using the locally supplemented Gridap.jl and GridapEmbedded.jl. Prior to following the below installation instructions ensure that no Manifest.toml file is present in the root folder of the project. In addition, ensure that either Julia 1.9 or 1.10 is used. Higher versions of Julia have not been tested.

1. Go to the root folder of the project and run the following command:
   ```
   julia --project=.
   ``` 
2. Then, enter the Pkg REPL by pressing `]`. Run the following command:
   ```
   dev ./deps/Gridap.jl ./deps/GridapEmbedded.jl
   ``` 
   The `Manifest.toml` file has now been generated in the root folder. Exit (first the PKG REPL and next) Julia by pressing backspace and typing `exit()`.
3. You can now proceed to run the test scripts by using the following command in the terminal (from the root folder)
   ```
   julia --project=. test/FIG8-1_h-convergence-on-a-cut-square.jl 
   ``` 
   where you can substitute the filename with any of the available test scripts. The results of the test will be printed to the screen. Note that computing the first result of each method ('std', 'BGP' or 'AL-BGP') takes some additional computational time (because of Julia's JIT compilation).

**IMPORTANT NOTE**
The results shown in the paper were generated with the PETSc interface of the [MUMPS](https://mumps-solver.org/index.php) sparse direct solver as provided by the Julia package [GridapPETSc](https://github.com/gridap/GridapPETSc.jl), as we observed it to be more numerically robust than the UMFPACK sparse direct solver bundled with Julia. In order to use GridapPETSc along with MUMPS, one needs to compile the PETSc library. The installation instructions of PETSc are machine-specific. Hence, we refer to the installation instructions of GridapPETSc and the associated github actions [CI.yml](https://github.com/gridap/GridapPETSc.jl/blob/master/.github/workflows/ci_extra.yml) file for additional details. Once GridapPETSc/PETSc/MUMPS has been properly installed on your system, then you need to set the `PETSc` variable in the test scripts to `true`. Otherwise, the UMFPACK sparse direct solver bundled with Julia will be used. 
