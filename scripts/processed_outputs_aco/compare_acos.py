import pandas as pd
import os

# Caminhos das pastas contendo os arquivos CSV
pasta_com_rvns = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/scripts/processed_outputs_aco/ACO-com-RVNS"
pasta_sem_rvns = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/scripts/processed_outputs_aco/ACO-sem-RVNS"

# Colunas de interesse
colunas = [
    "Graph Name", "$V(G)$", "$E(G)$", "$\delta$", "$\Delta$", 
    "Melhor Fitness", "Lower Bound", "Upper Bound", "Graph Density"
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
dados_com_rvns = carregar_dados(pasta_com_rvns)
dados_sem_rvns = carregar_dados(pasta_sem_rvns)

# Combinar os dados em uma única tabela
tabela_combinada = []
for graph_name in dados_com_rvns.keys():
    if graph_name in dados_sem_rvns:
        linha_com_rvns = dados_com_rvns[graph_name]
        linha_sem_rvns = dados_sem_rvns[graph_name]
        
        # Criar uma nova linha com os dados de ambos os testes
        nova_linha = {
            "Graph Name": graph_name,
            "$V(G)$": linha_com_rvns["$V(G)$"],
            "$E(G)$": linha_com_rvns["$E(G)$"],
            "Graph Density": linha_com_rvns["Graph Density"],
            "LB": linha_com_rvns["Lower Bound"],
            "UB": linha_com_rvns["Upper Bound"],
            "Menor Fitness (com RVNS)": linha_com_rvns["Melhor Fitness"],
            "Menor Fitness (sem RVNS)": linha_sem_rvns["Melhor Fitness"]        
        }
        tabela_combinada.append(nova_linha)

# Converter a lista combinada em um DataFrame
df_combinado = pd.DataFrame(tabela_combinada)

# Salvar a tabela combinada em um novo arquivo CSV
df_combinado.to_csv("tabela_combinada.csv", index=False)

print("Tabela combinada salva em 'tabela_combinada.csv'.")
