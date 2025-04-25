# Inexact Krylov methods for response calculations in density functional theory

[![][doi-badge]]([doi-link])
[![][arxiv-badge]]([arxiv-link])
[![][code-badge]]([code-link])

[doi-badge]: https://img.shields.io/badge/DOI-10.48550/arXiv.xxxx.xxxxx-blue
[doi-link]: https://doi.org/10.48550/arXiv.xxxx.xxxxx

[arxiv-badge]: https://img.shields.io/badge/arxiv-xxxx.xxxxx-red
[arxiv-link]: https://arxiv.org/abs/xxxx.xxxxx

[code-badge]: https://img.shields.io/badge/Julia-v1.10.3-blue.svg?logo=julia
[code-link]: https://julialang.org/downloads/

This repository contains the code to reproduce the results of the following paper:

Michael F. Herbst and Bonan Sun  
*Efficient Krylov methods for linear response in plane-wave electronic structure calculations*  
arXiv prerint, [arXiv:xxxx.xxxxx]([arxiv-link]), 2025, https://doi.org/10.48550/arXiv.xxxx.xxxxx

## Organization

Each directory contains the source code (and data, figures for Al40 and Fe2MnAl) for the corresponding system. For example, the structure of `Al40/` is as follows and other directories are similar.
```
Al40
 ┣ data_logs
 ┃ ┣ Al_extracted.jld2    # saved results
 ┃ ┗ summary_tables.md    # summary of the results for Al40
 ┣ figures                # figures
 ┃ ┗ repeat10             # figures for length 10 Al supercell
 ┃ ┃ ┣ Al.eps
 ┃ ┃ ┣ ...
 ┣ Al40.jl                # main source for Al systems
 ┣ GeneratePlots.jl       # generate figures for Al40
 ┣ GenerateTables.jl      # generate tables for Al40
 ┗ test.jl                # testing script for a small system
```

## Getting started

### Installing Julia

The recommended way to install Julia is to install [`juliaup`](https://github.com/JuliaLang/juliaup), a Julia version manager. See the [`juliaup`](https://github.com/JuliaLang/juliaup) repository or [Julia official website](https://julialang.org/install/) for more instructions.

Once you have `juliaup` installed, we recommend using the same julia version as the one used to develop this code (1.10.3) to avoid any potential issues. Simply open a terminal from the directory of this repository and run the following command to start Julia with the correct version and activate the project environment:
```bash
juliaup add 1.10.3
julia +1.10.3 --project
```

### Installing dependencies

Then run the following command in the Julia REPL activated above to install the dependencies:
```julia
using Pkg
Pkg.instantiate()
```
This will install all the required packages. Notably it will install **not** the main distribution of [DFTK](https://dftk.org), but [a branch of my fork]([dftk-reproduce-fork]). Though the methods proposed in the paper will be incorporated in the main distribution of DFTK in the future, this repository is meant to be self-contained and reproducible. This code may not work with the main distribution of DFTK.

[dftk-reproduce-fork]: https://github.com/bonans/DFTK.jl/tree/reproduce-paper

### Running the tests
Run a small test with
```julia
include("Al40/test.jl")
```
to see if everything is working. It should create a directory `repeat3/` under `Al40/data_logs/` with 3 `.log` files and 3 `.jld2` files. It should also create one plot in `Al40/figures/repeat3/`. It should finish in less than 10 minutes on a modern laptop. 

## Reproducing the results in the paper

### Generating the tables and figures in the paper
The tables and figures in the paper can be generated from saved results in `Al40/data_logs/Al_extracted.jld2` and `Fe2MnAl/data_logs/Fe2MnAl_extracted.jld2`. 
Run
```julia
include("Al40/GeneratePlots.jl"); include("Al40/GenerateTables.jl")
include("Fe2MnAl/GeneratePlots.jl"); include("Fe2MnAl/GenerateTables.jl")
```
to generate all the tables and figures in the paper (except for Table 4 for silicon) by loading the saved results `.jld2` files in `Al40/data_logs/` and `Fe2MnAl/data_logs/` since re-running all the tests would take thousands of core hours, see next section if you really want to do that or if you want to run part of the tests shown in the paper. See the table below for a description of the them.

| File name | Description |
|--------|-------------|
| **For Aluminum:** |      |
| [`summary_tables.md`](Al40/data_logs/summary_tables.md) | Tab. 2, SM1, SM2, SM3 in paper and supplement: summary of the results for Al systems of different sizes |
| [`convergence_naive_20_-9.svg`](Al40/figures/repeat10/convergence_naive_20_-9.svg) | Fig. 2 in paper: Convergence of the baseline methods for Al40 |
| [`convergence_20_-9.svg`](Al40/figures/repeat10/convergence_20_-9.svg) | Fig. 3 in paper: Convergence of our methods and baseline methods for Al40 |
| [`CG_convergence_20_-9.svg`](Al40/figures/repeat10/CG_convergence_20_-9.svg) | Fig. 3 in paper: Convergence of our methods and baseline methods in terms of the no. of Hamiltonian applications for Al40 |
| [`Al.eps`](Al40/figures/repeat10/Al.eps) | Fig. 3 in paper: a visualization of the Al40 system |
| [`CG_niters_20_-9.svg`](Al40/figures/repeat10/CG_niters_20_-9.svg) | Fig. SM3 in supplement: evolution of no. of CG iterations and tolerances during GMRES for Al40 |
| [`normdeltaV.svg`](Al40/figures/repeat10/normdeltaV.svg) | Fig. SM2 in supplement:  Hxc term values for Al40 |
|  **For Fe2MnAl:**  |     |
| [`summary_tables.md`](Fe2MnAl/data_logs/summary_tables.md) | Tab. 3 in paper: summary of the results for Fe2MnAl |
| [`Hconvergence_10_-9.svg`](Fe2MnAl/figures/Hconvergence_10_-9.svg) | Fig. 4 in paper: analogous to [`convergence_20_-9`](Al40/figures/repeat10/convergence_20_-9.svg) for Fe2MnAl |
| [`HCG_convergence_10_-9.svg`](Fe2MnAl/figures/HCG_convergence_10_-9.svg) | Not in paper: analogous to [`CG_convergence_20_-9`](Al40/figures/repeat10/CG_convergence_20_-9.svg) for Fe2MnAl |
| [`HCG_niters_10_-9.svg`](Fe2MnAl/figures/HCG_niters_10_-9.svg) | Not in paper: analogous to [`CG_niters_20_-9.svg`](Al40/figures/repeat10/CG_niters_20_-9.svg) for Fe2MnAl |
| [`CG_err_res.svg`](Fe2MnAl/figures/CG_err_res.svg) | Fig. SM1 in supplement: CG residuals and errors for the worst conditioned Sternheimer equation for Al40 and Fe2MnAl |

For the Silicon system, run
```julia
include("silicon/test.jl"); include("silicon/GetTable.jl")
```
to reproduce the Tab. 4 in the paper from scratch. It will take 2~3 hours on a modern laptop.

### Reproducing the results

To rerun the calculations for the metallic systems, let's take the Heusler alloy as an example. First run 
```julia
include("Fe2MnAl/Fe2MnAl.jl")
``` 
to import procedures and functions. Then run 
```julia
setup_model(debug)
```
to run the self-consistent field (SCF) calculation to set up the model, either in debugging mode (`debug="debug"`) or in production mode (`debug="full"`), where the former uses a small cutoff and the latter uses the same parameters as in the paper. After the SCF, the response calculation can be run with
```julia
run_gmres(; debug, restart, tol, adaptive, CG_tol_scale_choice, precon)
```
where:
- `debug`: boolean, default `false`. If `true`, run the calculation for the small system otherwise for the same system as in the paper.
- `restart`: integer, default `20`. The restart parameter for GMRES.
- `tol`: float, default `1e-12`. The tolerance for GMRES convergence.
- `adaptive`: boolean or string, default `true`. Possible values are `true` (adaptive strategies proposed in the paper) or `"D10"`, `"D100"`, `"D10_n"` (baseline methods, see Table 1 in the paper for the definition).
- `CG_tol_scale_choice`: string, default `"agr"`. The strategy to choose the CG tolerances in case of `adaptive=true`. Possible values are `"agr"`, `"hdmd"`, `"grt"` (see Table 1 for the definition, `"hdmd"` refers to the `bal` strategy in the paper).
- `precon`: boolean, default `false`. If `true`, use the Kerker preconditioner.

Note that one GMRES run will generate a `.log` file and a `.jld2` file in the `data_logs/` directory. The `.log` file contains the output of the GMRES run, while the `.jld2` file contains the results. Each would take **5~30 hours** in `debug="full"` mode depending on the choice of the above parameters, e.g., how tight the tolerance is, which strategy is used etc. 

For the aluminum system, similar procedures can be followed with the exception that one more parameter `repeat` is needed in both `setup_model` and `run_gmres` to specify the size of the supercell.

## Contact
This is research code and therefore not necessarily user-friendly and actively maintained. If you have questions or encounter problems, get in touch!


