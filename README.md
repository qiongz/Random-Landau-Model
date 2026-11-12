# Random Landau Model

Chern numbers of the eigenstates of two models of the integer quantum Hall plateau transition: the continuum Landau-level model and the disordered Hofstadter lattice model. In both, a short-range disorder potential is projected into the lowest Landau level or the lowest magnetic subband. The same pipeline imposes generalized boundary conditions, diagonalizes the projected Hamiltonian, and computes the Chern number of every eigenstate.

A nonzero Chern number marks a conducting state. The Chern numbers in one sample sum to $+1$, the Chern number of that level. Their count $N_c$ is the quantity scaled in Zhu, Wu, Bhatt, and Wan, Phys. Rev. B **99**, 024205 (2019), [arXiv:1804.00398](https://arxiv.org/abs/1804.00398).

The previous programs are kept on the `legacy` branch.

## Build

Requires C++17, OpenMP, and oneMKL (LP64, GNU OpenMP threading). The same link line is used on Intel and on AMD EPYC. `icpc` is not required.

```bash
cmake -S . -B build -DMKLROOT=/opt/intel/oneapi/mkl/2025.2
cmake --build build -j
./build/rlm check
```

`check` is the numerical test a new build or a new eigensolver has to pass:

- Continuum model, $N_\phi=16$, $16N_\phi$ scatterers, boundary-condition grid $L_g=16$, seed 1: the relative residual at one boundary condition is below $10^{-8}$, and the Chern numbers of the sample sum to $1$ within $10^{-3}$
- Lattice model, $3\times 3$, flux $\phi_0/3$ per plaquette, lowest subband, same grid and seed: the same two tests

The GPU eigensolver is optional. On this machine the working build is HIP (`gfx1201`, rocSOLVER). The CUDA path is the same source file and needs `nvcc`.

```bash
cmake -S . -B build -DMKLROOT=/opt/intel/oneapi/mkl/2025.2 -DRLM_GPU=HIP
cmake --build build -j
./build/rlm check
```

With a GPU solver compiled in, `check` also diagonalizes one projected matrix with both solvers and requires $\max|\Delta E|<10^{-8}$.

The models and the Chern-number convention are in [docs/physics.md](docs/physics.md). Flags, output columns, and timing are in [docs/usage.md](docs/usage.md).

## Parameters

`--mesh`, `--samples`, `--seed`, `--solver`, `--threads`, and `--out` are shared. `--mesh` is the paper's grid size $L_g$.

| Flag | Model | Meaning | Default |
| --- | --- | --- | --- |
| `--nphi` | continuum | Flux quanta $N_\phi$ in the sample | 16 |
| `--quanta` | continuum | Factor in the sample side $\sqrt{N_\phi\times\mathrm{quanta}}$ | 1 |
| `--impurity` | continuum | Scatterers per flux quantum. The paper uses 16 | 16 |
| `--off-head` | continuum | Momentum cutoff. Negative selects $\lfloor\sqrt{N_\mathrm{imp}}/2\rfloor+1$ | $-1$ |
| `--lx`, `--ly` | lattice | Lattice dimensions. `lx` must be divisible by `q` | 3, 3 |
| `--q` | lattice | Flux per plaquette is $\phi_0/q$. The paper uses $q=3$ | 3 |
| `--band` | lattice | Magnetic subband. `0` is the lowest | 0 |
| `--mesh` | both | Boundary-condition grid $L_g$ | 16 |
| `--samples` | both | Disorder realizations | 1 |
| `--seed` | both | Seed of `std::mt19937` | 1 |
| `--solver` | both | `mkl` or `gpu` | `mkl` |
| `--threads` | both | OpenMP threads over $\theta_2$. `0` uses the OpenMP default | 0 |

## Commands

```bash
# Continuum Landau-level model. The paper uses 16 scatterers per flux quantum
# and Lg = 50.
./build/rlm run continuum --nphi 16 --mesh 50 --samples 2 --seed 1 --threads 16

# Hofstadter lattice, flux phi_0/3, disorder projected into the lowest subband.
./build/rlm run lattice --lx 12 --ly 12 --q 3 --band 0 --mesh 30 --seed 1

# Same projected matrix, GPU eigensolver. The binary must be built with -DRLM_GPU.
./build/rlm run continuum --nphi 64 --mesh 16 --seed 1 --solver gpu
```

Each sample prints its total Chern number on standard error. Energies and Chern numbers go to `--out`, which defaults to `energy_Chern_No.txt`.

## Layout

```text
include/rlm/rlm.hpp    parameters, models, eigensolver
src/continuum.cpp      lowest-Landau-level projection
src/lattice.cpp        lowest-subband projection of the Hofstadter model
src/pipeline.cpp       sweep of the boundary conditions
src/chern.cpp          Fukui-Hatsugai-Suzuki Chern numbers
src/solver_mkl.cpp     LAPACKE_zheevd
src/solver_gpu.cu      rocSOLVER or cuSOLVER zheevd
src/cli.cpp            command line
src/check.cpp          regression
```

The projected Hamiltonian is assembled on the host from one formula. `--solver` replaces only the eigensolver. MKL sweeps $\theta_2$ with OpenMP, and each `zheevd` runs single-threaded. The GPU solver handles those diagonalizations serially and copies the eigenvectors back before the Chern numbers are formed.

## License

[MIT](https://choosealicense.com/licenses/mit/)
