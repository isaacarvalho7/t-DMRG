# Adiabatic_quench_hdf5_salve

As implementações feitas nesse ramo foram construídas usando a biblioteca ITensor na versão C++V3 [https://itensor.org/]. Os códigos foram construídos com o objetivo de avaliar o comportamento da função de onda de modelos 1D durante o quench de tempo finito.

Breve descrição do código:

*Linha 27--72*
Input para tamanho da cadeia, número de ocupação, valor das interações iniciais e finais, número de sweeps, pinning field, valor da parâmetro de relaxamento tau_U e configurações para evolução pós quench de varredura.

*Linha 79--183*
rotinas para configurar o preenchimento, o hamiltoniano do modelo, os valores de U(t) e V(t) do quench de varredura, e número mínimo de pontos para o gráfico. 

*Linha 186--205*
Configuração do que será medido no quench de varredura.

*Linha 209--1127*
Rotina para o quench de varredura com as medias do observações: entropia vN, correlações de carga e spin, valores locais para densidade de carga e de spin. 

*Linha 1134--1135*
Rotina para salvar função de onda ao término do quench de varredura e o valor dos sites. 

Neste código não está incluso as rotinas de evolução e medida pós quench de varredura. 






