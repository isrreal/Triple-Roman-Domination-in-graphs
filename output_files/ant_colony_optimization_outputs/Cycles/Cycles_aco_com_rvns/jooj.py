import pandas as pd
import os
import numpy as np

# Diretório contendo os arquivos CSV
diretorio_csv = "/Cycles_aco_com_rvns"  # Substitua pelo caminho correto

# Lista para armazenar os dados processados
dados_processados = []

number_of_ants = 1
iterations = 5
# Processar cada arquivo CSV no diretório
for arquivo in os.listdir(diretorio_csv):
    if arquivo.endswith(".csv"):
        caminho_arquivo = os.path.join(diretorio_csv, arquivo)
        df = pd.read_csv(caminho_arquivo)

        # Calcular as métricas para cada arquivo
        menor_fitness = df['fitness_1_5'].min()
        media_fitness = df['fitness_1_5'].mean()
        desvio_padrao_fitness = df['fitness_1_5'].std()
        tempo_menor_fitness = df.loc[df['fitness_1_5'].idxmin(), 'elapsed time(seconds)']
        media_tempo = df['elapsed time(seconds)'].mean()
        desvio_padrao_tempo = df['elapsed time(seconds)'].std()

        # Adicionar os resultados à lista
        dados_processados.append({
                'Graph Name': df['graph_name'].unique()[0],
                '$V(G)$': df['graph_order'].unique()[0],
                '$E(G)$': df['graph_size'].unique()[0],
                '$\delta$': df['graph_min_degree'].unique()[0],
                '$\Delta$': df['graph_max_degree'].unique()[0],
                'Menor Fitness': df[fitness_col].min(),
                'Média Fitness': df[fitness_col].mean(),
                'Desvio Padrão Fitness': df[fitness_col].std(),
                'Tempo Menor Fitness': df.loc[df[fitness_col].idxmin(), time_col],
                'Média Tempo': df[time_col].mean(),
                'Desvio Padrão Tempo': df[time_col].std(),
                'Lower Bound': df['lower_bound'].unique()[0],
                'Upper Bound': df['upper_bound'].unique()[0],
                'Graph Density': df['graph_density'].unique()[0],
            })

# Criar um DataFrame com os dados processados
df_resultados = pd.DataFrame(dados_processados)

# Salvar o resultado em um novo arquivo CSV
df_resultados.to_csv("melhores_iteracoes.csv", index=False)
print("Arquivo 'melhores_iteracoes.csv' gerado com sucesso!")
