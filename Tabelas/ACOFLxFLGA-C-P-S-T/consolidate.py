import pandas as pd
import os

# Caminhos das pastas contendo os arquivos CSV
pasta_aco = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/ACOFLxFLGA-C-P-S-T/Stars/ACO"
pasta_ga = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/ACOFLxFLGA-C-P-S-T/Stars/GA"

# Colunas de interesse
colunas = [
    "Graph Name", "$V(G)$", "$E(G)$", "Graph Density", "LB", "UB", "Melhor Fitness", "Tempo do Melhor Fitness"
]

# Função para carregar os dados de uma pasta
def carregar_dados(pasta):
    dados = {}
    for arquivo in os.listdir(pasta):
        if arquivo.endswith(".csv"):
            caminho = os.path.join(pasta, arquivo)
            df = pd.read_csv(caminho)
            for _, row in df.iterrows():
                graph_name = row["Graph Name"]
                dados[graph_name] = row[colunas]
    return dados

# Carregar dados das duas pastas
dados_aco = carregar_dados(pasta_aco)
dados_ga = carregar_dados(pasta_ga)

# Combinar os dados em uma única tabela
tabela_combinada = []
for graph_name in dados_aco.keys():
    if graph_name in dados_ga:
        linha_aco = dados_aco[graph_name]
        linha_ga = dados_ga[graph_name]
        
        # Criar uma nova linha com os dados de ambos os algoritmos
        nova_linha = {
            "Graph Name": graph_name,
            "$V(G)$": linha_aco["$V(G)$"],
            "$E(G)$": linha_aco["$E(G)$"],
            "Graph Density": linha_aco["Graph Density"],
            "LB": linha_aco["LB"],
            "UB": linha_aco["UB"],
            "Melhor Fitness (ACO-FL)": linha_aco["Melhor Fitness"],
            "Tempo do Melhor Fitness (ACO-FL)": linha_aco["Tempo do Melhor Fitness"],
            "Melhor Fitness (FL-GA)": linha_ga["Melhor Fitness"],
            "Tempo do Melhor Fitness (FL-GA)": linha_ga["Tempo do Melhor Fitness"]
        }
            
        tabela_combinada.append(nova_linha)

# Converter a lista combinada em um DataFrame
df_combinado = pd.DataFrame(tabela_combinada)

# Salvar a tabela combinada em um novo arquivo CSV
df_combinado.to_csv("comparacao_ACO_GA.csv", index=False)

print("Tabela combinada salva em 'comparacao_ACO_GA.csv'.")
