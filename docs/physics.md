# Models and Chern numbers

This program computes the eigenstates used to extract the localization length exponent of the integer quantum Hall plateau transition. The calculation follows Zhu, Wu, Bhatt, and Wan, Phys. Rev. B **99**, 024205 (2019), [arXiv:1804.00398](https://arxiv.org/abs/1804.00398). Two microscopic models are projected into one band of Chern number $`+1`$. A state with a nonzero Chern number is a conducting state. The number of conducting states per sample,

```math
N_c=\int_{-\infty}^{\infty}\rho_c(E)\,\mathrm{d}E,
```

scales as $`N_c\sim N_\phi^{1-1/(2\nu)}`$. The fit in that paper gives $`\nu=2.48\pm 0.02`$. This program writes the boundary-condition-averaged energy and the Chern number of every eigenstate. It does not count $`N_c`$ or fit $`\nu`$.

The transition is controlled by the localization length $`\xi`$, called the correlation length in Eq. (1) of the paper,

```math
\xi(E)\sim|E-E_c|^{-\nu}.
```

After the projection used here, particle-hole symmetry fixes the critical energy at $`E_c=0`$.

## Generalized boundary conditions

Both models are solved on an $`L_1\times L_2`$ sample with generalized boundary conditions

```math
t(L\hat{e}_j)\,|\psi\rangle=e^{i\theta_j}|\psi\rangle,
\qquad \theta_j\in[0,2\pi).
```

The boundary-condition space is the torus $`[0,2\pi)^2`$. The program divides each side into $`L_g`$ intervals (`--mesh`). The paper uses $`L_g=50`$ for every continuum size it reports, $`L_g=30`$ for lattices with $`L\le 30`$, and $`L_g\approx\sqrt{2/3}\,L`$ for $`L=36`$ to $`81`$. Energies stored by the program are averages over the grid. Chern numbers use the plaquettes of that grid, so only two adjacent rows of eigenvectors are kept.

## Continuum Landau-level model

The continuum Hamiltonian is a Landau level plus a short-range impurity potential. With the magnetic length $`l_B=1`$,

```math
N_\phi=\frac{L_1 L_2}{2\pi}
```

is the number of flux quanta. The paper takes square samples and $`N_\mathrm{imp}=16N_\phi`$ scatterers,

```math
V(\mathbf{r})=\sum_{i=1}^{N_\mathrm{imp}}\delta(\mathbf{r}-\mathbf{r}_i),
```

and projects $`V`$ into the lowest Landau level. In the stripe basis of that level the projected matrix element is

```math
\begin{aligned}
H_{jk}(\theta_1,\theta_2)
&=\sum_{m,n}V_{mn}\,
e^{-\pi(m^2+n^2)/(2N_\phi)}\,
e^{-i\pi mn/N_\phi}\\
&\quad\times
e^{i(m\theta_2-n\theta_1)/N_\phi}\,
e^{-i 2\pi m j/N_\phi}\,
\langle\phi_j|\phi_{k-n}\rangle.
\end{aligned}
```

On a square sample the Gaussian depends on the aspect ratio only through $`L_2/L_1`$, so the factor printed above is the one in Sec. II.2 of the paper. The overlap is nonzero only when the orbital index differs from $`k`$ by the Fourier index $`n`$, modulo $`N_\phi`$.

The program keeps a finite window of Fourier components, of half-width

```math
m_0=\lfloor\sqrt{N_\mathrm{imp}}/2\rfloor+1
```

unless `--off-head` is set. $`N_\mathrm{imp}`$ is $`\lfloor\texttt{impurity}\times N_\phi/\texttt{quanta}\rfloor`$, and the default concentration is the paper's value 16. Each scatterer is given an independent sign $`\pm 1`$. The paper writes the $`\delta`$-functions with equal weight; the signs are what the legacy CPU program used. Positions are drawn uniformly in a square of side $`\sqrt{N_\phi\times\texttt{quanta}}`$. The default `--quanta 1` gives side $`\sqrt{N_\phi}`$. The Fourier phase evaluated by the program is $`2\pi(m-m_0)x+2\pi(n-m_0)y`$ with those coordinates, again the legacy CPU phase. The paper's wavevector on a sample of side $`L=\sqrt{2\pi N_\phi}`$ is $`2\pi m/L`$.

The random numbers are `std::mt19937` with seed `--seed`. Sample $`k`$ first discards the draws that samples $`0,\ldots,k-1`$ would have used, so it does not depend on which earlier samples were diagonalized. This stream is not the legacy MKL VSL stream, so energies do not reproduce a legacy output file bit for bit. The sum of the Chern numbers does.

## Hofstadter lattice model

The lattice Hamiltonian is a square-lattice tight-binding model in a uniform magnetic field, plus on-site disorder,

```math
\mathcal{H}
=-\sum_{\langle i,j\rangle}\!\left(t\,e^{i\theta_{ij}}c_i^\dagger c_j+\mathrm{h.c.}\right)
+\sum_i\epsilon_i\,c_i^\dagger c_i.
```

The flux per plaquette is $`\phi/\phi_0=1/q`$. For the paper's choice $`q=3`$ the clean spectrum has three magnetic subbands, with Chern numbers $`+1`$, $`-2`$, and $`+1`$, and

```math
N_\phi=L_1 L_2/3
```

states in each. The program requires $`L_1`$ divisible by $`q`$ and projects the on-site potential into subband `--band`. Band 0 is the lowest. The clean subband orbitals are the eigenvectors of the $`q\times q`$ Hofstadter matrix at the twisted crystal momentum; the projected matrix is

```math
H_{ab}
=\sum_{\mathbf{r}}
\psi_a^*(\mathbf{r})\,\epsilon(\mathbf{r})\,\psi_b(\mathbf{r}).
```

The disorder variables $`\epsilon(\mathbf{r})`$ are drawn uniformly from $`[-1,1]`$. That is the paper's choice $`W=1`$. The projected calculation sets the hopping to $`t=0`$, also the paper's choice: once $`t/W\le 0.2`$ the localization-length exponent no longer depends on $`t`$, and $`t=0`$ makes the clean subband flat, so particle-hole symmetry pins $`E_c`$ at zero. Eigenvectors are rotated back to the lattice before the Chern numbers are computed, because the link variables are overlaps of those lattice wave functions.

On the lattice, $`\theta=2\pi`$ is the same boundary condition as $`\theta=0`$. The program diagonalizes $`\theta_2=0,\ldots,L_g-1`$ and copies the $`\theta_2=0`$ eigenvectors and eigenvalues onto the periodic image. The continuum model evaluates the endpoint $`\theta=2\pi`$ as its own boundary condition.

## Chern number

The boundary-condition average of the Hall conductance of eigenstate $`m`$ is the Chern number,

```math
\langle\sigma_{xy}(m)\rangle=C(m)\,\frac{e^2}{h}.
```

The program uses the link variables of Fukui, Hatsugai, and Suzuki. The Berry phase on one plaquette is the argument of the product of the four overlaps around it. Dividing each overlap by its modulus does not change that argument, so the program skips the per-link normalization. The argument is taken in $`(-\pi,\pi]`$, which is the branch whose failures the paper analyzes in Sec. III.3.

The two models then apply a different constant so that the printed numbers sum to $`+1`$:

- Continuum: add $`1/N_\phi`$ to every state.
- Lattice: change the sign of the plaquette sum, and do not add $`1/N_\phi`$.

These constants are properties of the two legacy implementations. The printed column is the Chern number whose sample sum the paper requires to be $`+1`$.

If the grid is coarse, a plaquette phase can fall outside $`(-\pi,\pi]`$ and the sample sum jumps by an integer. For $`N_\phi=16`$ and $`N_\phi=64`$ with $`L_g=16`$, the sum differs from 1 by about $`10^{-15}`$. For $`N_\phi=256`$ with $`L_g=4`$, the sum is 2. The paper rejects an entire sample whose total Chern number is not $`+1`$. This program keeps the sample and prints the sum on standard error; the run is not usable for $`N_c`$ until that sum is 1.

A conducting state has a nonzero Chern number. Localized states contribute 0. Because the sample sum is $`+1`$, the usual pattern is a single state with $`C=+1`$, but a perturbation can create a pair with $`+1`$ and $`-1`$ away from the band center. Those pairs barely change $`N_c`$ and move the width of $`\rho_c(E)`$, which is why the paper treats $`N_c`$ as the more stable observable. The output is not rounded to integers.

## Eigensolver

`diagonalize` has the same contract as `LAPACKE_zheevd`: a column-major Hermitian matrix with both triangles filled, `uplo = U`, eigenvalues in ascending order, and eigenvector $`j`$ in column $`j`$. MKL and the GPU solver consume that same host matrix. The GPU path is `rocsolver_zheevd` under HIP or `cusolverDnZheevd` under CUDA. Chern numbers are always formed from the eigenvectors after they are back on the host.
