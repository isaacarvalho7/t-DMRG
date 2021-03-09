# t_DMRG

As implementações feitas nesse ramo fora construídas usando a biblioteca ITensor na versão C++V3 [https://itensor.org/]. Os códigos foram construídos com o objetivo de avaliar o comportamento da função de onda de modelos 1D após um quench. Uma breve descrição das implementações é apresentada a seguir.

<Spinless.cc>

Contém a implemetação do hamiltoniano do modelo Spinless que avalia a densidade de carga local antes e após o quench. Está classe apresenta os operadores de campo, criação de aniquilação, e o operador de densidade local.

<hubbard.cc>

Contém uma implementação do modelo de hubbard estendido 1D que avalia, após o quench, o comportamento do parâmetro de ordem CDW no tempo, descrito na referência [http://dx.doi.org/10.1016/j.physb.2017.09.001]. Além disso, essa implementação oferece oferece a possibilidade de recomeçar o cálculo a partir a da última evolução, visto que os a função de onda é salva a cada intervalo de tempo previamente definido.  








