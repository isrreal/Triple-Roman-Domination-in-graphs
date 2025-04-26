import pandas as pd
import os

# Diretório contendo os arquivos CSV
input_dir = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/output_files/outputs/genetic_algorithm_outputs/general-outputs/Miscellaneous-Networks/mixedHeuristics/"  # Substitua pelo caminho real
output_file = os.path.join(input_dir, "consolidated_best_results.csv")

# Lista para armazenar os melhores resultados
best_results = []

# Percorre todos os arquivos CSV na pasta
for arquivo in os.listdir(input_dir):
    if arquivo.endswith(".csv"):
        caminho_arquivo = os.path.join(input_dir, arquivo)

        if os.path.getsize(caminho_arquivo) > 0:
            try:
                df = pd.read_csv(caminho_arquivo)

                # Coluna de fitness (identificar automaticamente)
                fitness_cols = [col for col in df.columns if "fitness_heuristic" in col]
                if not fitness_cols:
                    print(f"Erro: Nenhuma coluna de fitness encontrada em {arquivo}")
                    continue

                fitness_col = fitness_cols[0]  # Assume apenas uma heurística por arquivo
                time_col = "elapsed_time(seconds)"

                # Obtém a melhor linha (menor fitness)
                best_row = df.loc[df[fitness_col].idxmin()]

                # Coletando os dados
                graph_data = {
                    "Graph Name": best_row["graph_name"],
                    "$V(G)$": best_row["graph_order"],
                    "$E(G)$": best_row["graph_size"],
                    "Graph Density": round(best_row["graph_density"], 4),
                    "LB": best_row["lower_bound"],
                    "UB": best_row["upper_bound"],
                    "Melhor Fitness": best_row[fitness_col],
                    "Tempo do Melhor Fitness": round(best_row[time_col], 4)
                }

                best_results.append(graph_data)

            except Exception as e:
                print(f"Erro ao processar {arquivo}: {e}")

# Salvar os resultados em CSV
if best_results:
    df_final = pd.DataFrame(best_results).round(3)
    df_final.to_csv(output_file, index=False)
    print(f"Arquivo consolidado salvo em: {output_file}")
else:
    print("Nenhum dado foi processado.")

