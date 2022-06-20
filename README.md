# Adiabatic_quench

Descrever as funções do código (básicos) - Hubbard, cada repositório tem

As implementações feitas nesse ramo foram construídas usando a biblioteca ITensor na versão C++V3 [https://itensor.org/].



## Breve descrição do código sweep_quench.cc

Neste código estão definidas as rotinas para o quench de varredura e para evolução do sistema em um tempo posterior. É necessário configurar o arquivo inputfile para definir os parâmetros do cálculo, as quantidades medidas e se as funções de onda serão salvas. Deve-se escolher os parâmetros de interação para que ao menos |Uf-U0| ou |Vf-V0| seja diferente de zero. Neste ramo, além do código, estão inclusos os Makefile para compilação do código no cluster e na chacobo. 











