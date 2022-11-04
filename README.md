# Time-dependent DMRG (t-DMRG)

This code describes the one dimensional extend Hubbard model (EHM) at the equilibrium and non-equilibrium behaviors. 
We use Density Matrix Renormalization Group (DMRG)[1] within tensor networks [2] and its temporal extension, namely time-dependent DMRG (tDMRG). 
The code evaluates the correlations functions and local densities for spin e charger, entanglement entropy for bipartite systems, fidelity between ground and evolved state, energy quench and energy variance error. 

## Directories

1. Doc : Directory containing the Design Document 
2. Ed_hubbard: Directory containing the exact diagonalization of EHM. 
3. Main: Directory containing the main code sweep_quench.cc that is the main module of the t-DMRG code.


## Compilation

The libraries needed to compile the code are:

1. GNU compilers ( version >= 7.5.0)

2. OpenBLas (https://www.openblas.net/)

3. HDF5 (https://www.hdfgroup.org/solutions/hdf5/)

4. Itensor3.0 (https://itensor.org/)

## Documentation

The documentation of the main module (sweep_quenches.cc) as well as other implementation details can be obtained from the Design Document available in the Doc directory.

## Parallelization

The code is parallelized for shared memory architectures using OpenMP instructions (https://www.openmp.org/).
To control this level of optimization it is recommended to make use of the option:

export OMP_NUM_THREADS = number of threads


## Description of the time-dependent DMRG code

Neste código estão definidas as rotinas para o quench de varredura e para evolução do sistema em um tempo posterior. É necessário configurar o arquivo inputfile para definir os parâmetros do cálculo, as quantidades medidas e se as funções de onda serão salvas. Deve-se escolher os parâmetros de interação para que ao menos |Uf-U0| ou |Vf-V0| seja diferente de zero. Neste ramo, além do código, estão inclusos um exemplo do Makefile de compilação do código. 

## References

[1] SCHOLLWÖCK, Ulrich. The density-matrix renormalization group in the age of matrix product states. Annals of physics, v. 326, n. 1, p. 96-192, 2011.

[2] https://itensor.org/












