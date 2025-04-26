import pandas as pd
import os

# Conjunto de grafos
graph_set = [
    "BAI" #"Harwell-Boeing", "Miscellaneous-Networks", "random_graphs", "BAI", "Harwell-Boeing", "Miscellaneous-Networks", "Trees", "Cycles", "Paths", "Stars"
]

# Parâmetros
number_of_ants = 1
iterations = 5

# Diretório de saída
output_dir = "processed_outputs_aco"
os.makedirs(output_dir, exist_ok=True)

for current_graph_set in graph_set:
    caminho_pasta = f'/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/output_files/ant_colony_optimization_outputs/atmost333_outputs/{current_graph_set}'
    dados_processados = []
    
    for arquivo in os.listdir(caminho_pasta):
        if arquivo.endswith('.csv'):
            caminho_arquivo = os.path.join(caminho_pasta, arquivo)
            df = pd.read_csv(caminho_arquivo)

            # Processamento
            fitness_col = f'fitness_{number_of_ants}_{iterations}'
            time_col = 'elapsed_time(seconds)'
            
            entrada = {
                'Graph Name': df['graph_name'].unique()[0],
                '$V(G)$': df['graph_order'].unique()[0],
                '$E(G)$': df['graph_size'].unique()[0],
                'Graph Density': round(df['graph_density'].unique()[0], 4),
                'LB': df['lower_bound'].unique()[0],
                'UB': df['upper_bound'].unique()[0],
                'Média de Fitness': round(df[fitness_col].mean(), 4),
                'Desvio Padrão Fitness': round(df[fitness_col].std(), 4),
                'Média de Tempo': round(df[time_col].mean(), 4),
                'Desvio Padrão de Tempo': round(df[time_col].std(), 4)
            }
            
            dados_processados.append(entrada)
    
    # Salvar o DataFrame processado
    df_final = pd.DataFrame(dados_processados)
    save_path = os.path.join(output_dir, f'ACO_{current_graph_set}.csv')
    df_final.to_csv(save_path, index=False)
    print(f"Arquivo salvo: {save_path}")

