# Time-dependent DMRG (t-DMRG)

Codigo que ..... para o modelo de Hubbard extendido (EHM) ....
Descrever as funções do código (básicos) - Hubbard, cada repositório tem

## Diretorios

1. Doc : Diretorio contendo Design Document 
2. Ed_hubbard: Diretorio contendo codigo para diagonalizacao exata (ED) do modelo de Hubbard extendido.
3. Principal: Contem o codigo sweep_quench.cc que o modulo principal do codigo t-DMRG.


## Compilacao

As bibliotecas necessarias para a compilacao do codigo sao:

1. Compiladores GNU ( versao >= 7.5.0)

2. OpenBLas (https://www.openblas.net/)

3. HDF5 (https://www.hdfgroup.org/solutions/hdf5/)

4. Itensor3.0 (https://itensor.org/)

## Documentacao 

A documentacao do modulo principal (sweep_quenches.cc) bem como outros detalhes da implementacao podem ser obtidos no Design Documento disponivel no diretorio doc.

## Paralelizacao 

O codigo foi paralelizado para arquitetura de memoria compartilhada via instrucoes OpenMP (https://www.openmp.org/).
Para controlar esse nivel de otimizacao e recomendado fazer uso da opcao:

export OMP_NUM_THREADS = numero de threads


## Breve descrição do código time-dependent DMRG

Neste código estão definidas as rotinas para o quench de varredura e para evolução do sistema em um tempo posterior. É necessário configurar o arquivo inputfile para definir os parâmetros do cálculo, as quantidades medidas e se as funções de onda serão salvas. Deve-se escolher os parâmetros de interação para que ao menos |Uf-U0| ou |Vf-V0| seja diferente de zero. Neste ramo, além do código, estão inclusos os Makefile para compilação do código no cluster e na chacobo. 


## Referencias 

[1] Paper principal DMRG

[2] Citacao para lib Itensor












