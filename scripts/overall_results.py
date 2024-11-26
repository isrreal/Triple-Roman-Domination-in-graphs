import pandas as pd
import os
import sys
# Caminho da pasta com os arquivos CSV

pasta = sys.argv[1]
caminho_pasta = f'../output_files/genetic_algorithm_tests/{pasta}'

# Lista para armazenar os dados processados
dados_processados = []

# Itera sobre cada arquivo CSV na pasta
for arquivo in os.listdir(caminho_pasta):
    if arquivo.endswith('.csv'):
        caminho_arquivo = os.path.join(caminho_pasta, arquivo)
        df = pd.read_csv(caminho_arquivo)

        # Aplica as operações nas colunas desejadas
        entrada = {
            'Graph Name': df['graph_name'].unique()[0],
            '$V(G)$': df['graph_order'].unique()[0],
            '$E(G)$': df['graph_size'].unique()[0],
            '$\delta$': df['graph_min_degree'].unique()[0],
            '$\Delta$': df['graph_max_degree'].unique()[0],
            'Lower Bound': df['lower_bound'].unique()[0],
            'Upper Bound': df['upper_bound'].unique()[0],
            'Fitness Value': df['GA_fitness_heuristic2'].min(),
            'Elapsed Time(seconds)': df['elapsed_time_GA(seconds)'].min(),
        }
        
        dados_processados.append(entrada)

# Converte para um DataFrame e salva em um único CSV
df_final = pd.DataFrame(dados_processados)
df_final.to_csv('saida_unificada.csv', index=False)




#import pandas as pd
#import sys
#import matplotlib as plt
#
#csv_file = sys.argv[1]
#dataframe = pd.read_csv(csv_file)
#
##labels = dataframe[0:]
#
## Extrai informações e corrige o comprimento das listas
#graph_name = dataframe['graph_name']
#
##print(labels)
#
#data = {
#    'Graph name': graph_name, 
#    '$V(G)$': dataframe['graph_order'],
#    '$E(G)$': dataframe['graph_size'],
#    '$\delta(G)$': dataframe['graph_min_degree'], 
#    '$\Delta(G)$': dataframe['graph_max_degree'],
#   'Lower Bound': dataframe['lower_bound'],
#    'Upper Bound': dataframe['upper_bound'],
#   'Minimum Fitness Value': dataframe['GA_fitness_heuristic1'].min(),
#  'Minimum Elapsed Time': dataframe['elapsed_time_GA(seconds)'].min()
#}
#
#dataframe_table = pd.DataFrame(data)
#
#latex_table = dataframe_table.to_latex(index=False,
#    caption="Comparação do Algoritmo Genético sem RVNS no conjunto representativo", 
#    position="htbp",  # The preferred positions where the table should be placed in the document ('here', 'top', 'bottom', 'page')
#    column_format="|c|c|c|c|c|c|c|c|c|",  # The format of the columns: left-aligend first column and center-aligned remaining columns as per APA guidelines
#    escape=False,  # Disable escaping LaTeX special characters in the DataFrame
#    float_format="{:0.2f}".format  # Formats floats to two decimal places)
#)
#
#print(latex_table)
#
