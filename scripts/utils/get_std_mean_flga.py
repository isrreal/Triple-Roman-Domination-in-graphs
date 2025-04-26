import pandas as pd
import os

# Conjunto de grafos e heurísticas
graph_set = [
	"random_graphs"#"BAI", "Harwell-Boeing", "Miscellaneous-Networks"
]

heuristics = {
    "4": "mixedHeuristics"
}

# Diretório de saída
output_dir = "processed_outputs"
os.makedirs(output_dir, exist_ok=True)

# Processamento
for current_graph_set in graph_set:
    # Criar uma subpasta para o conjunto de grafos atual
    graph_set_dir = os.path.join(output_dir, current_graph_set)
    os.makedirs(graph_set_dir, exist_ok=True)

    # Pasta principal dos arquivos CSV
    base_path = f"/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/output_files/outputs/genetic_algorithm_outputs/general-outputs/{current_graph_set}"

    # Verifica se o caminho base existe
    if not os.path.exists(base_path):
        continue

    # Dicionário para armazenar os dados consolidados por grafo
    graph_data = {}

    for heuristic, heuristic_folder in heuristics.items():
        heuristic_path = os.path.join(base_path, heuristic_folder)

        if not os.path.exists(heuristic_path):
            print(f"Pasta não encontrada: {heuristic_path}")
            continue

        for arquivo in os.listdir(heuristic_path):
            if arquivo.endswith('.csv'):
                caminho_arquivo = os.path.join(heuristic_path, arquivo)

                if os.path.getsize(caminho_arquivo) > 0:
                    try:
                        df = pd.read_csv(caminho_arquivo)
                        fitness_col = f'fitness_heuristic_{heuristic}'
                        time_col = 'elapsed_time(seconds)'
                        graph_name = df['graph_name'].unique()[0]

                        if graph_name not in graph_data:
                            graph_data[graph_name] = {
                                'Graph Name': graph_name,
                                '$V(G)$': df['graph_order'].unique()[0],
                                '$E(G)$': df['graph_size'].unique()[0],
                                'Graph Density': round(df['graph_density'].unique()[0], 4),
                                'LB': df['lower_bound'].unique()[0],
                                'UB': df['upper_bound'].unique()[0],
                                'Fitness Values': [],  # Lista para armazenar valores de fitness
                                'Time Values': []     # Lista para armazenar valores de tempo
                            }

                        # Adicionar valores de fitness e tempo ao dicionário
                        graph_data[graph_name]['Fitness Values'].extend(df[fitness_col].tolist())
                        graph_data[graph_name]['Time Values'].extend(df[time_col].tolist())

                    except Exception as e:
                        print(f"Erro ao processar {caminho_arquivo}: {e}")

    # Calcular médias e desvios padrão para cada grafo
    best_results = []
    for graph_name, data in graph_data.items():
        fitness_values = data['Fitness Values']
        time_values = data['Time Values']

        mean_fitness = round(sum(fitness_values) / len(fitness_values), 4) if fitness_values else None
        std_fitness = round(pd.Series(fitness_values).std(), 4) if fitness_values else None
        mean_time = round(sum(time_values) / len(time_values), 4) if time_values else None
        std_time = round(pd.Series(time_values).std(), 4) if time_values else None

        best_results.append({
            'Graph Name': graph_name,
            '$V(G)$': data['$V(G)$'],
            '$E(G)$': data['$E(G)$'],
            'Graph Density': data['Graph Density'],
            'LB': data['LB'],
            'UB': data['UB'],
            'Média de Fitness': mean_fitness,
            'Desvio Padrão Fitness': std_fitness,
            'Média de Tempo': mean_time,
            'Desvio Padrão de Tempo': std_time
        })

    # Salvar os resultados consolidados
    df_final = pd.DataFrame(best_results)
    save_path = os.path.join(graph_set_dir, f'GA_{current_graph_set}_best_fitness_comparison.csv')
    df_final.to_csv(save_path, index=False)
    print(f"Arquivo salvo: {save_path}")
