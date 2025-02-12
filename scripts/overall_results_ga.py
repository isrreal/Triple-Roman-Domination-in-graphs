import pandas as pd
import os
import sys

# Conjunto de grafos e heurísticas
graph_set = [
    "BAI", "Harwell-Boeing", "Miscellaneous-Networks", "random_graphs", "Trees", "Cycles", "Paths", "Stars"
]
heuristics = ["1", "2", "3", "4"]

# Diretório de saída
output_dir = "processed_outputs"
os.makedirs(output_dir, exist_ok=True)

# Processamento
dados_processados = []
for current_graph_set in graph_set:
    for current_heuristic in heuristics:
        caminho_pasta = f'../output_files/genetic_algorithm_outputs/{current_graph_set}/{current_heuristic}'
        if not os.path.exists(caminho_pasta):
            continue
        
        for arquivo in os.listdir(caminho_pasta):
            if arquivo.endswith('.csv'):
                caminho_arquivo = os.path.join(caminho_pasta, arquivo)
                df = pd.read_csv(caminho_arquivo)

                # Processamento
                fitness_col = f'fitness_heuristic_{current_heuristic}'
                time_col = 'elapsed_time(seconds)'
                
                entrada = {
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
                }
                
                dados_processados.append(entrada)
    
    # Salvar o DataFrame processado
    df_final = pd.DataFrame(dados_processados)
    save_path = os.path.join(output_dir, f'GA_{current_graph_set}.csv')
    df_final.to_csv(save_path, index=False)
    print(f"Arquivo salvo: {save_path}")

