# Usage

The models, the Chern number, and the sum rule are defined in [physics.md](physics.md). This page is the command line and the regression.

## Running

```bash
./build/rlm run continuum --nphi 32 --impurity 16 --quanta 1 --mesh 50 \
  --samples 4 --seed 1 --threads 16 --solver mkl --out energy_Chern_No.txt

./build/rlm run lattice --lx 12 --ly 12 --q 3 --band 0 --mesh 30 \
  --samples 4 --seed 1 --threads 16
```

`--help` prints the flags and returns 0. A missing or unknown flag returns 2. A failed diagonalization, or an output file that cannot be written, returns 1.

The continuum model rejects $N_\phi<1$, $L_g<1$, a non-positive `--quanta`, fewer than one scatterer, or a momentum cutoff that is not strictly smaller than $N_\phi$. The lattice model requires positive `lx`, `ly`, `q`, and `mesh`, `lx` divisible by `q`, and `band` in $[0,q)$.

## Output

`--out` has one row per eigenstate, in scientific notation with 16 significant digits:

```text
sample    energy    Chern number
```

The energy is the average over the $(L_g+1)^2$ boundary conditions that the sweep stores. Standard error prints one line per sample:

```text
continuum sample 0 chern_sum 1.0000000000000000e+00
```

If that sum differs from 1 by more than about $10^{-3}$, the grid is too coarse for this $N_\phi$. Increase `--mesh`. Do not count conducting states from that sample. The paper discards any sample whose Chern numbers do not add to $+1$.

Count a state as conducting when its Chern number is nonzero. Most such states are near $+1$. A $+1,-1$ pair away from $E=0$ is still two conducting states, and the sample sum remains $+1$. The program does not round or threshold the column.

## Threads and eigensolvers

On the MKL path, OpenMP distributes the $\theta_2$ points. Each `zheevd` sets the local MKL thread count to 1, so the library does not spawn more threads inside the OpenMP region. `--threads 0`, or omitting the flag, uses the OpenMP default. `--threads 1` walks $\theta_2$ serially.

`--solver gpu` exists only in a binary configured with `-DRLM_GPU=HIP` or `-DRLM_GPU=CUDA`. Those diagonalizations run serially on one device; `--threads` does not split them. Without the GPU solver compiled in, the flag reports that the binary must be rebuilt. The HIP architecture defaults to `gfx1201` and is changed with `-DRLM_HIP_ARCH`.

## check

```bash
./build/rlm check
```

`check` does not write `energy_Chern_No.txt`. Success ends with a line `PASS` and exit code 0. Any failed test is printed as `FAIL` and the exit code is 1.

| Model | Sample | Tests |
| --- | --- | --- |
| Continuum | $N_\phi=16$, 16 scatterers per flux quantum, $L_g=16$, seed 1, sample 0 | Relative residual at $\theta=(0,0)$ below $10^{-8}$; energies finite; Chern numbers sum to $1$ within $10^{-3}$ |
| Lattice | $3\times 3$, $q=3$, lowest subband, $L_g=16$, seed 1, sample 0 | The same three tests |

The residual is $\lVert H\psi-E\psi\rVert/\lVert E\psi\rVert$, and the reported value is the maximum over the eigenstates at that one boundary condition.

If a GPU solver was compiled in, the same projected matrix at $\theta=(0,0)$ is diagonalized by both solvers, for one continuum sample and one lattice sample, and $\max|\Delta E|$ must be below $10^{-8}$. A continuum sample with $L_g=4$ is also run on both solvers and is required only to return finite numbers. Per-state Chern numbers on that coarse grid are not a failure condition.

## Timing

One sample performs $(L_g+1)^2$ dense diagonalizations of an $N_\phi\times N_\phi$ Hermitian matrix. On this EPYC 9V74, with 16 OpenMP threads:

| Sample | MKL | HIP |
| --- | --- | --- |
| Continuum, $N_\phi=64$, $L_g=16$ | 0.12 s | 1.32 s |
| Lattice, $12\times 12$, $q=3$ ($N_\phi=48$), $L_g=16$ | 0.08 s | |

Both totals of the Chern number are 1. At this matrix size the transfer onto the GPU costs more than `zheevd`. After changing the eigensolver, compare eigenvalues in `check` before increasing $L_g$ until the Chern numbers of a sample sum to 1.
