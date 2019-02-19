# Random-Landau-Model
Diagonalization and Chern number calculation of each eigenstate for the disorder projected continuum Landau-level model and Hofstadter lattice model
## Description
With Chern number we could distinguish the conducting and total density of states and thus perform finite-size scaling and extract the localization-length exponent. This program calculate the disorder projected continuum random Landau-level model (rlm) as well as the Hofstadter lattice model (hlm), diagonalizing the projected hamiltonian and calculate the Chern No. for each eigenstates in each disorder sample. This program support both CPU/GPU or mixed parallel acceleration.For more information of the two models please visit APS: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.024205/ or arXiv:https://arxiv.org/abs/1804.00398
## Prerequisites
Environment:Linux

Libraries:
* icc,mkl:https://software.intel.com/en-us/parallel-studio-xe
* g++,g++-multilib,liblapack-dev,libblas-dev
* CUDA toolkit 9.0:https://developer.nvidia.com/cuda-90-download-archive?target_os=Linux
* MAGMA library:http://icl.cs.utk.edu/projectsfiles/magma/doxygen/installing.html

## Install and compile
ubuntu:

>sudo apt-get install g++ g++-multilib build-essential liblapack-dev libblas-dev

Intel icc,mkl are necessary to be installed:
>make cpu

To use GPU version, CUDA are necessary, and when compiling, the default g++ should be switched to g++6.0:
> make rlm_gpu_cusolver rlm_gpu_mkl

The most efficient CPU/GPU combined diagonalization library is the MAGMA library, with MAGMA installed:
> make rlm_gpu_magma
## Running and benchmark
### Hofstadter lattice model (hlm)
> make hlm_cpu_dp
<pre><code>./hlm_cpu_dp -h
./hlm_cpu_dp: option requires an argument -- 'h'
Usage: ./hlm_cpu_dp -[OPTIONS]
Options:
  -x  lx
  -y  ly
  -q  q [p/q, p set to 1]
  -g  theta mesh griding size (n_mesh)
  -b  band to be projected onto
  -c  # of disorder configurations(n_sample)
  -s  seed to generate seeds
  -t  # of threads
</code></pre>

For an example, here we diagonalize a 3x3 lattice model with 1/q (q=3) flux quantum per unit cell,and calculate Chern No. with grid size g=64, project the disorder potential to  the lowest subband of the three (b=0),and generate 2 disorder samples, using 4 threads:

> ./hlm_cpu_dp -x 3 -y 3 -q 3 -g 64 -b 0 -c 2 -t 4  

and the results are stored in "energy_Chern_No.txt" in three columns, the 1st column is the disorder sample #,the 2nd column is the twist boundary condition averaged eigenenergy, and the 3rd column is the corresponding Chern #:

<pre><code>more energy_Chern_No.txt
       0     -0.42254796 1.4918622e-16
       0    -0.071045388 1
       0       0.2389036 4.3042826e-17
       1     -0.76097491 -1.2045486e-16
       1     -0.52535084 2.1239521e-16
       1     -0.13183907 1 </code></pre>

As a rule of thumb, the sum of Chern # in each sample is "1".  
### Random Landau model (rlm)
#### CPU version
> make rlm_cpu_dp
<pre><code>./rlm_cpu_dp -h
./rlm_cpu_dp: option requires an argument -- 'h'
Usage: ./rlm_cpu_dp -[OPTIONS]
Options:
  -n  n_phi
  -q  quanta concentration
  -i  impurity concentration
  -g  theta mesh griding size (n_mesh)
  -d  disorder potential K-space cut (offhead)
  -c  # of disorder configurations(n_sample)
  -s  seed to generate seeds
  -t  # of threads </code></pre>

Similar to the hlm, here we diagonalize a system with total flux quanta n_phi=16, quanta concentration=1,impurity concentration=16, and mesh grid size g=64, generate 2 disordered samples and using 4 threads:  

>./rlm_cpu_dp -n 16 -i 16 -g 64 -c 2 -t 4  

And the results are written in "energy_Chern_No.txt", we could verify that there are two samples with each 16 eigenstate and each sample has total Chern # equals "1".    

<pre><code>more energy_Chern_No.txt  
0        -5.7376113      -1.2490009e-16
0        -4.7726244      6.1756156e-16
0        -4.5293731      -6.2450045e-16
0        -3.8107133      -4.0245585e-16
0        -3.1913483      -8.3266727e-17
0        -2.4917556      1
0        -1.8941804      3.8857806e-16
0        -1.3713557      6.0368377e-16
0        -0.63845787     3.469447e-17
0        0.071196911     -1.2490009e-16
0        0.54396831      -3.6082248e-16
0        1.4574762       2.220446e-16
0        3.4363649       -2.7755576e-17
0        4.9878217       1.3877788e-17
0        7.3543744       -2.7755576e-17
0        10.586218       -1.3877788e-17
1        -7.6545345      9.7144515e-17
1        -6.3143188      1.3877788e-17
1        -5.4981977      1.3877788e-17
1        -4.7521262      9.7144515e-17
1        -3.6033109      1.2490009e-16
1        -2.9968285      1.8041124e-16
1        -2.2325506      -7.6327833e-16
1        -1.7562629      1.7416624e-15
1        -1.2725596      -1.9012569e-15
1        -0.52287515     1
1        0.34612954      3.5388359e-16
1        1.0266966       -3.3306691e-16
1        2.0576746       -2.0816682e-16
1        3.0926          -2.7755576e-17
1        3.7359859       1.110223e-16
1        6.3444782       0  Landau-Model
</code></pre>    

#### GPU version and benchmark
The GPU version is originally used for single precision calculation but finally I found that the accuracy of the Chern # calculation is severely effected for nphi>1500, double precision CPU version is suggested for finite-size-scaling calculation of localization-length exponent, the single precision CPU/GPU version are still useful for generating and store wave function for further Machine learning test.

Here we test the CPU+GPU version with MAGMA library supported:
>make rlm_gpu_magma

benchmark it with nphi=3000, mesh grid size g=4,and # of threads 4:
<pre><code>time ./rlm_gpu_magma -n 3000 -g 4 -t 4
coefficients initialization time:       0.0322808 s
potential initialization time:          0.260278 s
diagonalization time/k-point:           1.12105 s
Chern No. calculation time/k-point:     0.0242928 s
estimated total time:                   000:00:28
estimated # of samples/hour:            124.98
# module                                 %(time)
coefficients                            0.112068
potential                               0.903599
diagonalization                         97.2976
Chern No. calculation                   1.68673

real	0m27.645s
user	1m59.882s
sys	0m7.497s
</code></pre>

we compare it with single precision CPU version:
> make rlm_cpu_sp
<pre><code>time ./rlm_cpu_sp -n 3000 -g 4 -t 4
coefficients initialization time:       0.004526 s
potential initialization time:          19.6461 s
diagonalization time/k-point:           4.05126 s
Chern No. calculation time/k-point:     0.0192458 s
estimated total time:                   000:02:01
estimated # of samples/hour:            29.6743
# module                                %(time)
coefficients                            0.00373072
potential                               16.194
diagonalization                         83.485
Chern No. calculation                   0.317281


real	2m4.944s
user	5m33.394s
sys	0m1.096s</code></pre>

Here the CPU+GPU version takes 27.645 seconds and CPU version takes 124.944 seconds, the speedup is 4.51x, with the hardware configuration gtx1060,i7-8650u. We could see 97% of the time in the GPU version is spent on diagonalization, since the potential initialiation function could be written in a CUDA kernel function and there's no bottleneck on this. When nphi is larger, the speedup could be much better.
## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
## License
[MIT](https://choosealicense.com/licenses/mit/)
