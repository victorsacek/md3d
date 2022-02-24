Sequência de operações para a criação do modelo teste levando-se em consideração variações composicionais:

>ipython gera_temper/gera_temper_3D.py
A partir do arquivo da espessura da litosfera Section A-A.txt, cria os seguintes arquivos 
(1) campo de temperatura inicial - Temper_0_3D.txt 
(2) interfaces composicionais - interfaces.txt
O arquivo interfaces.txt possui a seguinte estrutura:
1a linha: fator de escala da viscosidade para as n camadas do modelo
2a linha: densidade de cada camada: kg/m3
3a linha: produção de calor radiogênico de cada camada: W/kg
4a linha em diante: profundidade das (n-1) interfaces

>copiar os arquivos interfaces.txt e Temper_0_3D.txt para a pasta da execução do modelo.
(nessa versão do MD3D a leitura do arquivo de temperatura é feita diretamente do arquivo *.txt, não sendo necessário converter para binário.

>ajustar o arquivo param_1.5.3.txt (observar que o número de interface no arquivo param deve ser o mesmo do número de interfaces contidas no arquivo interfaces.txt

>ajustar o arquivo roda.sh para o número de processos necessários e não esquecer de incluir “-te 1”

>ipython tudo2D.py
Gera gráficos de Temperatura e outras propriedades do modelo ao longo do tempo.
