import pandas as pd
import os

# Caminhos das pastas contendo os arquivos CSV
pasta_aco = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/Tabelas Restantes/Tabela BAI/ACO/"
pasta_ga = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/Tabelas/Tabelas Restantes/Tabela BAI/GA/"
arquivo_pli = "/home/leggen-assis/Documents/faculdade/TCC/implementacao/3RDF/results/point-of-interest/resultados_BAI.csv"

# Colunas de interesse
colunas = [
    "Graph Name", "$V(G)$", "$E(G)$", "Graph Density", "UB", "Melhor Fitness", "Tempo do Melhor Fitness"
]

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

def carregar_pli(arquivo):
    df_pli = pd.read_csv(arquivo)
    dados_pli = {}
    for _, row in df_pli.iterrows():
        graph_name = row["filename"].replace(".txt", "")  # Remove a extensão .txt
        pli_value = row["TRD_number"]

        # Verifica se o valor é numérico e converte para inteiro se não houver casas decimais
        if pd.notna(pli_value):
            pli_value = int(pli_value) if pli_value == int(pli_value) else round(pli_value, 2)

        status = row["status"]
        dados_pli[graph_name] = f"{pli_value}*" if status == "Optimal" else pli_value
    return dados_pli

# Carregar dados das três fontes
dados_aco = carregar_dados(pasta_aco)
dados_ga = carregar_dados(pasta_ga)
dados_pli = carregar_pli(arquivo_pli)

# Combinar os dados em uma única tabela
tabela_combinada = []
for graph_name in dados_aco.keys():
    if graph_name in dados_ga:
        linha_aco = dados_aco[graph_name]
        linha_ga = dados_ga[graph_name]
        pli_value = dados_pli[graph_name]

        # Determinar o melhor fitness encontrado
        melhor_fitness = min(linha_aco["Melhor Fitness"], linha_ga["Melhor Fitness"])
        
        # Remover '*' e converter para float
        pli_value_numeric = float(str(pli_value).replace("*", ""))
        gap_relativo = abs(pli_value_numeric - melhor_fitness) * 100 / melhor_fitness

        nova_linha = {
            "Graph Name": graph_name,
            "$V(G)$": linha_aco["$V(G)$"],
            "$E(G)$": linha_aco["$E(G)$"],
            "Graph Density": round(linha_aco["Graph Density"], 3),
            "UB": linha_aco["UB"],
            "F/ACO-FL": linha_aco["Melhor Fitness"],
            "T/ACO-FL": round(linha_aco["Tempo do Melhor Fitness"], 3),
            "F/FL-GA": linha_ga["Melhor Fitness"],
            "T/FL-GA": round(linha_ga["Tempo do Melhor Fitness"], 3),
            "PLI": pli_value,  # Adicionada vírgula que estava faltando
            "Gap relativo": round(gap_relativo, 2)
        }
        
        tabela_combinada.append(nova_linha)

# Converter a lista combinada em um DataFrame
df_combinado = pd.DataFrame(tabela_combinada)

# Salvar a tabela combinada em um novo arquivo CSV
df_combinado.to_csv("comparacao_ACO_GA.csv", index=False)

print("Tabela combinada salva em 'comparacao_ACO_GA.csv'.")

